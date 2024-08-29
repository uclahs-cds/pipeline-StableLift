include { compress_and_index_HTSlib } from './utils.nf'

process predict_stability_StableLift {
    container params.docker_image_stablelift
    containerOptions "-v ${moduleDir}:${moduleDir}"

    publishDir path: "${params.output_dir_base}/output",
        pattern: "stability.tsv",
        mode: "copy",
        saveAs: { "StableLift-${sample_id}-${variant_caller}.tsv" }

    input:
    tuple val(sample_id), path(features_rds)
    path(rf_model)
    val(variant_caller)

    output:
    tuple val(sample_id), path("stability.tsv"), emit: stability_tsv

    script:
    """
    Rscript "${moduleDir}/scripts/predict-variant-stability.R" \
        --variant-caller "${variant_caller}" \
        --features-dt "${features_rds}" \
        --rf-model "${rf_model}" \
        --output-tsv "stability.tsv"
    """

    stub:
    """
    touch "stability.tsv"
    """
}

process run_apply_stability_annotations {
    container params.docker_image_bcftools

    publishDir path: "${params.output_dir_base}/output",
        pattern: "{stability,filtered}.vcf.gz{,.tbi}",
        mode: "copy",
        saveAs: { "${sample_id}-${it}" }

    input:
    tuple val(sample_id),
        path(annotated_vcf, stageAs: 'inputs/*'),
        // FIXME Should there be an annotated_vcf_tbi?
        path(stability_tsv, stageAs: 'inputs/*'),
        path(stability_tsv_tbi, stageAs: 'inputs/*')

    output:
    tuple val(sample_id),
        path("stability.vcf.gz"),
        path("stability.vcf.gz.tbi"),
        emit: stability_vcf_with_index
    tuple val(sample_id),
        path("filtered.vcf.gz"),
        path("filtered.vcf.gz.tbi"),
        emit: filtered_vcf_with_index

    script:
    stability_vcf = "stability.vcf.gz"
    stability_vcf_tbi = "${stability_vcf}.tbi"

    filtered_vcf = "filtered.vcf.gz"
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
    touch "stability.vcf.gz"
    touch "stability.vcf.gz.tbi"
    touch "filtered.vcf.gz"
    touch "filtered.vcf.gz.tbi"
    """
}

workflow workflow_predict_stability {
    take:
    vcf_with_sample_id
    r_annotations
    rf_model
    variant_caller

    main:

    predict_stability_StableLift(
        r_annotations,
        rf_model,
        variant_caller
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
