process run_extract_vcf_features {
    container params.docker_image_stablelift
    containerOptions "-v ${moduleDir}:${moduleDir}"

    input:
    tuple val(sample_id), path(vcf)

    output:
    tuple val(sample_id), path('output.Rds'), emit: r_annotations

    script:
    intermediate_filepath = "${params.output_dir_base}/stablelift-${params.stablelift_version}/intermediate/${task.process}"

    slug = "stablelift-${sample_id}"

    """
    Rscript "${moduleDir}/scripts/extract-vcf-features.R" \
        --input-vcf "${vcf}" \
        --variant-caller ${params.variant_caller} \
        --output-rds "output.Rds"
    """
}

workflow extract_features {
    take:
    vcf_with_sample_id

    main:
    if (params.variant_caller == "HaplotypeCaller") {
        error "HaplotypeCaller is not supported yet"
    } else {
        run_extract_vcf_features(vcf_with_sample_id)
    }

     emit:
     r_annotations = run_extract_vcf_features.out.r_annotations
}
