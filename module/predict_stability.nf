include { compress_and_index_HTSlib } from './utils.nf'

process predict_stability_StableLift {
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

    stub:
    """
    touch "${output_file_name}"
    """
}

process run_apply_stability_annotations {
    container params.docker_image_bcftools

    publishDir path: "${params.output_dir_base}/output",
        pattern: "*.vcf.gz{,.tbi}",
        mode: "copy"

    input:
    tuple val(sample_id),
        path(annotated_vcf, stageAs: 'inputs/*'),
        // FIXME Should there be an annotated_vcf_tbi?
        path(stability_tsv, stageAs: 'inputs/*'),
        path(stability_tsv_tbi, stageAs: 'inputs/*')

    output:
    tuple val(sample_id),
        path(stability_vcf),
        path(stability_vcf_tbi),
        emit: stability_vcf_with_index
    tuple val(sample_id),
        path(filtered_vcf),
        path(filtered_vcf_tbi),
        emit: filtered_vcf_with_index

    script:
    slug = "${sample_id}_LiftOver"

    stability_vcf = "${slug}_stability.vcf.gz"
    stability_vcf_tbi = "${stability_vcf}.tbi"

    filtered_vcf = "${slug}_filtered.vcf.gz"
    filtered_vcf_tbi = "${filtered_vcf}.tbi"

    """
    bcftools annotate \
        -a "${stability_tsv}" \
        -c CHROM,POS,STABILITY_SCORE,STABILITY \
        -h <(echo '##INFO=<ID=STABILITY_SCORE,Number=1,Type=String,Description="Proportion of trees in random forest model predicting `class = concordant`">
##INFO=<ID=STABILITY,Number=1,Type=String,Description="Stability status: STABLE or UNSTABLE">') \
        -o "${stability_vcf}" \
        "${annotated_vcf}"

    bcftools index -t "${stability_vcf}"

    bcftools filter \
        -i 'INFO/STABILITY="STABLE"' \
        -o "${filtered_vcf}" \
        "${stability_vcf}"

    bcftools index -t "${filtered_vcf}"
    """

    stub:
    """
    touch "${stability_vcf}"
    touch "${stability_vcf_tbi}"
    touch "${filtered_vcf}"
    touch "${filtered_vcf_tbi}"
    """
}

workflow workflow_predict_stability {
    take:
    vcf_with_sample_id
    r_annotations
    rf_model

    main:

    predict_stability_StableLift(
        r_annotations,
        rf_model
    )

    compress_and_index_HTSlib(
        predict_stability_StableLift.out.stability_tsv
    )

    run_apply_stability_annotations(
        vcf_with_sample_id.join(
            compress_and_index_HTSlib.out.compressed_tsv_with_index,
            failOnDuplicate: true,
            failOnMismatch: true
        )
    )

    emit:
    stability_vcf = run_apply_stability_annotations.out.stability_vcf_with_index
    filtered_vcf =  run_apply_stability_annotations.out.filtered_vcf_with_index
}
