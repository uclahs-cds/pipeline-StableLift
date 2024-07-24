
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

    slug = "${sample_id}"

    """
    bgzip ${tsv} --output output.tsv.gz

    tabix \
        --sequence 1 \
        --begin 2 \
        --end 2 \
        output.tsv.gz
    """

    stub:
    """
    touch "output.tsv.gz"
    touch "output.tsv.gz.tbi"
    """
}
