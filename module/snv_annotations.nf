include { compress_and_index_HTSlib } from './utils.nf'

process run_Funcotator_GATK {
    container params.docker_image_gatk

    publishDir path: "${params.output_dir_base}/GATK-${params.gatk_version}/intermediate/${task.process}",
        pattern: "output.vcf.gz",
        mode: "copy",
        enabled: params.save_intermediate_files,
        saveAs: { "Funcotator-${sample_id}-${dest_fasta_id}.vcf.gz" }

    input:
        tuple val(sample_id),
            path(vcf, stageAs: 'inputs/*'),
            path(index, stageAs: 'inputs/*')
        tuple val(dest_fasta_id), path(dest_fasta_ref), path(dest_fasta_fai), path(dest_fasta_dict)
        path (funcotator_sources)

    output:
        tuple val(sample_id), path('output.vcf.gz'), emit: funcotator_vcf

    script:
        """
        gatk Funcotator \
            --variant "${vcf}" \
            --reference "${dest_fasta_ref}" \
            --ref-version "${dest_fasta_id}" \
            --data-sources-path "${funcotator_sources}" \
            --output-file-format VCF \
            --output "output.vcf.gz"
        """

    stub:
    """
    touch "output.vcf.gz"
    """
}

process annotate_RepeatMasker_BCFtools {
    container params.docker_image_bcftools

    publishDir path: "${params.output_dir_base}/BCFtools-${params.bcftools_version}/intermediate/${task.process}",
        pattern: "output.vcf.gz",
        mode: "copy",
        enabled: params.save_intermediate_files,
        saveAs: { "RepeatMasker-${sample_id}.vcf.gz" }

    input:
    tuple val(sample_id), path(vcf, stageAs: 'inputs/*')
    path(repeat_bed)

    output:
    tuple val(sample_id), path('output.vcf.gz'), emit: repeatmasker_vcf

    script:

    """
    bcftools annotate \
        -a "${repeat_bed}" \
        -c CHROM,FROM,TO,REPEAT \
        -h <(echo '##INFO=<ID=REPEAT,Number=1,Type=String,Description="RepeatMasker Region">') \
        -o "output.vcf.gz" \
        "${vcf}"
    """

    stub:
    """
    touch "output.vcf.gz"
    """
}

process extract_TrinucleotideContext_BEDTools {
    container params.docker_image_bedtools

    publishDir path: "${params.output_dir_base}/BEDtools-${params.bedtools_version}/intermediate/${task.process}",
        pattern: "output.{bed,tsv}",
        mode: "copy",
        enabled: params.save_intermediate_files,
        saveAs: { "Trinucleotide-${sample_id}-${dest_fasta_id}.${file(it).getExtension()}" }

    input:
    tuple val(sample_id), path(vcf)
    tuple val(dest_fasta_id), path(dest_fasta_ref), path(dest_fasta_fai), path(dest_fasta_dict)

    output:
    tuple val(sample_id), path('output.tsv'), emit: trinucleotide_tsv
    tuple val(sample_id), path('output.bed'), emit: trinucleotide_bed

    script:
    """
    zcat ${vcf} \
        | grep -v "^#" \
        | awk '{if (length(\$4) == 1 && length(\$5) == 1) print \$1"\t"(\$2-2)"\t"(\$2+1)"\t"\$2}' \
        > "output.bed"

    bedtools getfasta \
        -fi ${dest_fasta_ref} \
        -bed "output.bed" \
        -fo "partial.tsv" \
        -name \
        -tab

    paste -d '\t' \
        <(cut -f 1,4 "output.bed") \
        <(cut -f 2 "partial.tsv") \
        > "output.tsv"
    """

    stub:
    """
    touch "output.tsv"
    touch "output.bed"
    """
}

process annotate_trinucleotide_BCFtools {
    container params.docker_image_bcftools

    publishDir path: "${params.output_dir_base}/BCFtools-${params.bcftools_version}/intermediate/${task.process}",
        pattern: "output.vcf.gz",
        mode: "copy",
        enabled: params.save_intermediate_files,
        saveAs: { "Trinucleotide-annotated-${sample_id}.vcf.gz" }

    input:
    tuple val(sample_id),
        path(vcf, stageAs: 'inputs/*'),
        path(tsv, stageAs: 'inputs/*'),
        path(tsv_tbi, stageAs: 'inputs/*')

    output:
    tuple val(sample_id), path('output.vcf.gz'), emit: trinucleotide_vcf

    script:
    """
    bcftools annotate \
        --annotations ${tsv} \
        --columns CHROM,POS,TRINUCLEOTIDE \
        --header-lines <(echo '##INFO=<ID=TRINUCLEOTIDE,Number=1,Type=String,Description="Trinucleotide Context">') \
        --output output.vcf.gz \
        ${vcf}
    """

    stub:
    """
    touch "output.vcf.gz"
    """
}

workflow workflow_apply_snv_annotations {
    take:
    vcf_with_sample_id
    dest_fasta_data

    main:

    run_Funcotator_GATK(
        vcf_with_sample_id,
        dest_fasta_data,
        Channel.value(params.funcotator_data.data_source)
    )

    annotate_RepeatMasker_BCFtools(
        run_Funcotator_GATK.out.funcotator_vcf,
        Channel.value(params.repeat_bed)
    )

    extract_TrinucleotideContext_BEDTools(
        annotate_RepeatMasker_BCFtools.out.repeatmasker_vcf,
        dest_fasta_data
    )

    compress_and_index_HTSlib(
        extract_TrinucleotideContext_BEDTools.out.trinucleotide_tsv
    )

    annotate_trinucleotide_BCFtools(
        annotate_RepeatMasker_BCFtools.out.repeatmasker_vcf.join(
            compress_and_index_HTSlib.out.compressed_tsv_with_index,
            failOnDuplicate: true,
            failOnMismatch: true
        )
    )

    emit:
    annotated_vcf = annotate_trinucleotide_BCFtools.out.trinucleotide_vcf
}
