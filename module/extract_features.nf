process run_extract_vcf_features {
    container params.docker_image_stablelift
    containerOptions "-v ${moduleDir}:${moduleDir}"

    publishDir path: "${intermediate_filepath}",
        pattern: "features.Rds",
        mode: "copy",
        enabled: params.save_intermediate_files,
        saveAs: { "${slug}.${file(it).getExtension()}" }

    input:
    tuple val(sample_id), path(vcf)

    output:
    tuple val(sample_id), path('features.Rds'), emit: r_annotations

    script:
    intermediate_filepath = "${params.output_dir_base}/stablelift-${params.stablelift_version}/intermediate/${task.process}"

    slug = "stablelift-${sample_id}"

    """
    Rscript "${moduleDir}/scripts/extract-vcf-features.R" \
        --input-vcf "${vcf}" \
        --variant-caller ${params.variant_caller} \
        --output-rds "features.Rds"
    """
}

process predict_variant_stability {
    container params.docker_image_stablelift
    containerOptions "-v ${moduleDir}:${moduleDir}"

    publishDir path: "${params.output_dir_base}/output",
        pattern: "${output_file_name}",
        mode: "copy"

    input:
    tuple val(sample_id), path(features_rds)
    path(rf_model)

    output:
    tuple val(sample_id), path(output_file_name), emit: stability_tsv

    script:
    output_file_name = "stablelift-${sample_id}.tsv"

    """
    Rscript "${moduleDir}/scripts/predict-liftover-stability.R" \
        --features-dt "${features_rds}" \
        --rf-model "${rf_model}" \
        --variant-caller "${params.variant_caller}" \
        --output-tsv "${output_file_name}"
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
        run_extract_vcf_features.out.r_annotations.set { ch_annotations }
    }

    predict_variant_stability(
        ch_annotations,
        Channel.value(params.rf_model)
    )

    emit:
    r_annotations = predict_variant_stability.out.stability_tsv
}
