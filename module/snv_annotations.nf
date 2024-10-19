include { compress_and_index_vcf; compress_and_index_tsv } from './utils.nf'

process add_genotype_field {
    container params.docker_image_stablelift
    containerOptions "-v ${moduleDir}:${moduleDir}"

    publishDir path: "${params.output_dir_base}/intermediate/${task.process.replace(':', '/')}",
        pattern: "output.vcf",
        mode: "copy",
        enabled: params.save_intermediate_files,
        saveAs: { "${sample_id}_add-GT.vcf" }

    input:
        tuple val(sample_id), path(vcf), path(index)

    output:
        tuple val(sample_id), path('output.vcf'), emit: vcf_with_gt

    script:
    """
    Rscript "${moduleDir}/scripts/add-GT-Strelka2.R" \
        --input-vcf "${vcf}" \
        --output-vcf output.vcf.gz
    gunzip output.vcf.gz
    """

    stub:
    """
    touch output.vcf
    """
}

process run_Funcotator_GATK {
    container params.docker_image_gatk

    publishDir path: "${params.output_dir_base}/intermediate/${task.process.replace(':', '/')}",
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
        def ref_version_map = [GRCh37: 'hg19', GRCh38: 'hg38']

        """
        gatk Funcotator \
            --variant "${vcf}" \
            --reference "${dest_fasta_ref}" \
            --ref-version "${ref_version_map[dest_fasta_id]}" \
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

    publishDir path: "${params.output_dir_base}/intermediate/${task.process.replace(':', '/')}",
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

    publishDir path: "${params.output_dir_base}/intermediate/${task.process.replace(':', '/')}",
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

    publishDir path: "${params.output_dir_base}/intermediate/${task.process.replace(':', '/')}",
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
    tuple val(sample_id), path('output.vcf.gz'), path('output.vcf.gz.tbi'), emit: trinucleotide_vcf

    script:
    """
    bcftools annotate \
        --annotations ${tsv} \
        --columns CHROM,POS,TRINUCLEOTIDE \
        --header-lines <(echo '##INFO=<ID=TRINUCLEOTIDE,Number=1,Type=String,Description="Trinucleotide Context">') \
        --write-index=tbi \
        --output output.vcf.gz \
        ${vcf}
    """

    stub:
    """
    touch "output.vcf.gz"
    touch "output.vcf.gz.tbi"
    """
}

workflow workflow_annotate_snvs {
    take:
    vcf_with_sample_id
    dest_fasta_data

    main:

    if (params.variant_caller == "Strelka2") {
        add_genotype_field(
            vcf_with_sample_id
        )
        compress_and_index_vcf(
            add_genotype_field.out.vcf_with_gt
        )
        vcf_for_annotation = compress_and_index_vcf.out.compressed_vcf_with_index
    } else {
        vcf_for_annotation = vcf_with_sample_id
    }

    run_Funcotator_GATK(
        vcf_for_annotation,
        dest_fasta_data,
        Channel.value(params.funcotator_data_source)
    )

    annotate_RepeatMasker_BCFtools(
        run_Funcotator_GATK.out.funcotator_vcf,
        Channel.value(params.repeat_bed)
    )

    extract_TrinucleotideContext_BEDTools(
        annotate_RepeatMasker_BCFtools.out.repeatmasker_vcf,
        dest_fasta_data
    )

    compress_and_index_tsv(
        extract_TrinucleotideContext_BEDTools.out.trinucleotide_tsv
    )

    annotate_trinucleotide_BCFtools(
        annotate_RepeatMasker_BCFtools.out.repeatmasker_vcf.join(
            compress_and_index_tsv.out.compressed_tsv_with_index,
            failOnDuplicate: true,
            failOnMismatch: true
        )
    )

    emit:
    annotated_vcf = annotate_trinucleotide_BCFtools.out.trinucleotide_vcf
}
