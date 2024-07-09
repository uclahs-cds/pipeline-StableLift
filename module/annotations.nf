process run_repeatmasker {
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
    intermediate_filepath = "${params.output_dir_base}/bcftools-${params.bcftools_version}/intermediate/${task.process}"

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

process run_trinucleotide_context {
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
    intermediate_filepath = "${params.output_dir_base}/bedtools-${params.bedtools_version}/intermediate/${task.process}"

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

process run_compress_and_index {
    container params.docker_image_samtools

    publishDir path: "${intermediate_filepath}",
        pattern: "output.tsv.gz",
        mode: "copy",
        enabled: params.save_intermediate_files

    publishDir path: "${intermediate_filepath}",
        pattern: "output.tsv.gz.tbi",
        mode: "copy",
        enabled: params.save_intermediate_files

    input:
    tuple val(sample_id), path(tsv)

    output:
    tuple val(sample_id), path('output.tsv.gz'), path('output.tsv.gz.tbi'), emit: compressed_tsv_with_index

    script:
    intermediate_filepath = "${params.output_dir_base}/samtools-${params.bcftools_version}/intermediate/${task.process}"

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


process run_trinucleotide_annotate {
    container params.docker_image_bcftools

    publishDir path: "${intermediate_filepath}",
        pattern: "output.vcf.gz",
        mode: "copy",
        enabled: params.save_intermediate_files,
        saveAs: { "${slug}.vcf.gz" }

    input:
    tuple
        val(sample_id),
        path(vcf, stageAs: 'inputs/*'),
        path(tsv, stageAs: 'inputs/*'),
        path(tsv_tbi, stageAs: 'inputs/*')

    // FIXME Should this process also emit the index file? It seems like
    // bcftools won't produce it without the --write-index flag
    output:
    tuple val(sample_id), path('output.vcf.gz'), emit: trinucleotide_vcf

    script:
    intermediate_filepath = "${params.output_dir_base}/bcftools-${params.bcftools_version}/intermediate/${task.process}"

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



workflow apply_annotations {
    take:
    vcf_with_sample_id
    dest_fasta_data

    main:
    run_repeatmasker(
        vcf_with_sample_id,
        Channel.value(params.repeat_bed)
    )

    run_trinucleotide_context(
        run_repeatmasker.out.repeatmasker_vcf,
        dest_fasta_data
    )

    run_compress_and_index(run_trinucleotide_context.out.trinucleotide_tsv)

     run_trinucleotide_annotate(
         run_repeatmasker.out.repeatmasker_vcf.join(
             run_compress_and_index.out.compressed_tsv_with_index,
             failOnDuplicate: true,
             failOnMismatch: true
         )
     )

     emit:
     annotated_vcf = run_trinucleotide_annotate.out.trinucleotide_vcf
}
