include { compress_and_index_tsv } from './utils.nf'

process predict_stability_StableLift {
    container params.docker_image_stablelift
    containerOptions "-v ${moduleDir}:${moduleDir}"

    input:
    tuple val(sample_id), path(features_rds)
    path(rf_model)
    val(variant_caller)

    output:
    tuple val(sample_id), path("stability.tsv"), emit: stability_tsv

    script:
    spec_arg = (params.getOrDefault('target_specificity', null) != null) ? "--specificity \"${params.get('target_specificity')}\"" : ""
    thresh_arg = (params.getOrDefault('target_threshold', null) != null) ? "--threshold \"${params.get('target_threshold')}\"" : ""

    """
    Rscript "${moduleDir}/scripts/predict-variant-stability.R" \
        --variant-caller "${variant_caller}" \
        --features-dt "${features_rds}" \
        --rf-model "${rf_model}" \
        --output-tsv "stability.tsv" \
        ${spec_arg} \
        ${thresh_arg}
    """

    stub:
    """
    touch "stability.tsv"
    """
}

process run_apply_stability_annotations {
    container params.docker_image_bcftools

    publishDir path: "${params.output_dir_base}/output",
        pattern: "StableLift-${dest_fasta_id}{,_filtered}.vcf.gz{,.tbi}",
        mode: "copy",
        saveAs: { "${sample_id}_${variant_caller}_${it}" }

    input:
    tuple val(sample_id),
        path(annotated_vcf, stageAs: 'inputs/*'),
        path(stability_tsv, stageAs: 'inputs/*'),
        path(stability_tsv_tbi, stageAs: 'inputs/*')
    val(variant_caller)
    val(dest_fasta_id)

    output:
    tuple val(sample_id),
        path("StableLift-${dest_fasta_id}.vcf.gz"),
        path("StableLift-${dest_fasta_id}.vcf.gz.tbi"),
        emit: stability_vcf_with_index
    tuple val(sample_id),
        path("StableLift-${dest_fasta_id}_filtered.vcf.gz"),
        path("StableLift-${dest_fasta_id}_filtered.vcf.gz.tbi"),
        emit: filtered_vcf_with_index

    script:
    stability_vcf = "StableLift-${dest_fasta_id}.vcf.gz"
    stability_vcf_tbi = "${stability_vcf}.tbi"

    filtered_vcf = "StableLift-${dest_fasta_id}_filtered.vcf.gz"
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
    touch "StableLift-${dest_fasta_id}.vcf.gz"
    touch "StableLift-${dest_fasta_id}.vcf.gz.tbi"
    touch "StableLift-${dest_fasta_id}_filtered.vcf.gz"
    touch "StableLift-${dest_fasta_id}_filtered.vcf.gz.tbi"
    """
}

workflow workflow_predict_stability {
    take:
    vcf_with_sample_id
    r_annotations
    rf_model
    variant_caller
    dest_fasta_id

    main:

    predict_stability_StableLift(
        r_annotations,
        rf_model,
        variant_caller
    )

    compress_and_index_tsv(
        predict_stability_StableLift.out.stability_tsv
    )

    run_apply_stability_annotations(
        vcf_with_sample_id.join(
            compress_and_index_tsv.out.compressed_tsv_with_index
            // failOnDuplicate: true,
            // failOnMismatch: true
        ),
        variant_caller,
        dest_fasta_id
    )

    emit:
    stability_vcf = run_apply_stability_annotations.out.stability_vcf_with_index
    filtered_vcf =  run_apply_stability_annotations.out.filtered_vcf_with_index
}
