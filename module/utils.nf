
process compress_and_index_HTSlib {
    container params.docker_image_samtools

    publishDir path: "${params.output_dir_base}/intermediate/${task.process.replace(':', '/')}",
        pattern: "output.tsv.gz{,.tbi}",
        mode: "copy",
        enabled: params.save_intermediate_files,
        saveAs: { "${sample_id}${file(it).getName().replace(file(it).getSimpleName(), "")}" }

    input:
    tuple val(sample_id), path(tsv)

    output:
    tuple val(sample_id),
        path('output.tsv.gz'),
        path('output.tsv.gz.tbi'),
        emit: compressed_tsv_with_index

    script:
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
