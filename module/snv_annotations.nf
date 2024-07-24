process run_Funcotator_GATK {
    container params.docker_image_gatk

    publishDir path: "${intermediate_filepath}",
        pattern: "output.vcf.gz",
        mode: "copy",
        enabled: params.save_intermediate_files,
        saveAs: { "${slug}.vcf.gz" }

    input:
        tuple val(sample_id),
            path(vcf, stageAs: 'inputs/*'),
            path(index, stageAs: 'inputs/*')
        tuple val(dest_fasta_id), path(dest_fasta_ref), path(dest_fasta_fai), path(dest_fasta_dict)
        path (funcotator_sources)

    output:
        tuple val(sample_id), path('output.vcf.gz'), emit: funcotator_vcf

    script:
        intermediate_filepath = "${params.output_dir_base}/GATK-${params.gatk_version}/intermediate/${task.process}"

        slug = "Funcotator-${sample_id}-${dest_fasta_id}"

        """
        gatk Funcotator \
            --variant "${vcf}" \
            --reference "${dest_fasta_ref}" \
            --ref-version "${dest_fasta_id}" \
            --data-sources-path "${funcotator_sources}" \
            --output-file-format VCF \
            --output "output.vcf.gz"
        """
}

process annotate_RepeatMasker_BCFtools {
    container params.docker_image_bcftools

    publishDir path: "${intermediate_filepath}",
        pattern: "output.vcf.gz",
        mode: "copy",
        enabled: params.save_intermediate_files,
        saveAs: { "${slug}.vcf.gz" }

    input:
    tuple val(sample_id), path(vcf, stageAs: 'inputs/*')
    path(repeat_bed)

    output:
    tuple val(sample_id), path('output.vcf.gz'), emit: repeatmasker_vcf

    script:
    intermediate_filepath = "${params.output_dir_base}/BCFtools-${params.bcftools_version}/intermediate/${task.process}"

    slug = "RepeatMasker-${sample_id}"

    """
    bcftools annotate \
        -a "${repeat_bed}" \
        -c CHROM,FROM,TO,REPEAT \
        -h <(echo '##INFO=<ID=REPEAT,Number=1,Type=String,Description="RepeatMasker Region">') \
        -o "output.vcf.gz" \
        "${vcf}"
    """
}

process extract_TrinucleotideContext_BEDTools {
    container params.docker_image_bedtools

    publishDir path: "${intermediate_filepath}",
        pattern: "output.bed",
        mode: "copy",
        enabled: params.save_intermediate_files,
        saveAs: { "${slug}.bed" }

    publishDir path: "${intermediate_filepath}",
        pattern: "output.tsv",
        mode: "copy",
        enabled: params.save_intermediate_files,
        saveAs: { "${slug}-full.tsv" }

    input:
    tuple val(sample_id), path(vcf)
    tuple val(dest_fasta_id), path(dest_fasta_ref), path(dest_fasta_fai), path(dest_fasta_dict)

    output:
    tuple val(sample_id), path('output.tsv'), emit: trinucleotide_tsv
    tuple val(sample_id), path('output.bed'), emit: trinucleotide_bed

    script:
    intermediate_filepath = "${params.output_dir_base}/BEDtools-${params.bedtools_version}/intermediate/${task.process}"

    slug = "Trinucleotide-${sample_id}"

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
}

process compress_and_index_HTSlib {
    container params.docker_image_samtools

    publishDir path: "${intermediate_filepath}",
        pattern: "output.tsv.gz{,.tbi}",
        mode: "copy",
        enabled: params.save_intermediate_files

    input:
    tuple val(sample_id), path(tsv)

    output:
    tuple val(sample_id), path('output.tsv.gz'), path('output.tsv.gz.tbi'), emit: compressed_tsv_with_index

    script:
    intermediate_filepath = "${params.output_dir_base}/SAMtools-${params.samtools_version}/intermediate/${task.process}"

    slug = "Trinucleotide-${sample_id}"

    """
    bgzip ${tsv} --output output.tsv.gz

    tabix \
        --sequence 1 \
        --begin 2 \
        --end 2 \
        output.tsv.gz
    """
}


process annotate_trinucleotide_BCFtools {
    container params.docker_image_bcftools

    publishDir path: "${intermediate_filepath}",
        pattern: "output.vcf.gz",
        mode: "copy",
        enabled: params.save_intermediate_files,
        saveAs: { "${slug}.vcf.gz" }

    input:
    tuple val(sample_id),
        path(vcf, stageAs: 'inputs/*'),
        path(tsv, stageAs: 'inputs/*'),
        path(tsv_tbi, stageAs: 'inputs/*')

    output:
    tuple val(sample_id), path('output.vcf.gz'), emit: trinucleotide_vcf

    script:
    intermediate_filepath = "${params.output_dir_base}/BCFtools-${params.bcftools_version}/intermediate/${task.process}"

    slug = "Trinucleotide-${sample_id}"

    """
    bcftools annotate \
        --annotations ${tsv} \
        --columns CHROM,POS,TRINUCLEOTIDE \
        --header-lines <(echo '##INFO=<ID=TRINUCLEOTIDE,Number=1,Type=String,Description="Trinucleotide Context">') \
        --output output.vcf.gz \
        ${vcf}
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
