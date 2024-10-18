process compress_and_index_vcf {
    container params.docker_image_samtools

    publishDir path: "${params.output_dir_base}/intermediate/${task.process.replace(':', '/')}",
        pattern: "output.vcf.gz{,.tbi}",
        mode: "copy",
        enabled: params.save_intermediate_files,
        saveAs: { "${sample_id}${file(it).getName().replace(file(it).getSimpleName(), "")}" }

    input:
    tuple val(sample_id), path(input_vcf)

    output:
    tuple val(sample_id),
        path('output.vcf.gz'),
        path('output.vcf.gz.tbi'),
        emit: compressed_vcf_with_index

    script:
    """
    bgzip ${input_vcf} --output output.vcf.gz
    tabix -p vcf output.vcf.gz
    """

    stub:
    """
    touch output.vcf.gz
    touch output.vcf.gz.tbi
    """
}

process compress_and_index_tsv {
    container params.docker_image_samtools

    publishDir path: "${params.output_dir_base}/intermediate/${task.process.replace(':', '/')}",
        pattern: "output.tsv.gz{,.tbi}",
        mode: "copy",
        enabled: params.save_intermediate_files,
        saveAs: { "${sample_id}${file(it).getName().replace(file(it).getSimpleName(), "")}" }

    input:
    tuple val(sample_id), path(input_tsv)

    output:
    tuple val(sample_id),
        path('output.tsv.gz'),
        path('output.tsv.gz.tbi'),
        emit: compressed_tsv_with_index

    script:
    """
    bgzip ${input_tsv} --output output.tsv.gz
    tabix -s 1 -b 2 -e 2 output.tsv.gz
    """

    stub:
    """
    touch output.tsv.gz
    touch output.tsv.gz.tbi
    """
}