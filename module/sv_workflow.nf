
process liftover_SV_StableLift{
    container params.docker_image_stablelift

    publishDir path: "${params.output_dir_base}/intermediate/${task.process.replace(':', '/')}",
        pattern: "liftover.vcf.gz",
        mode: "copy",
        enabled: params.save_intermediate_files,
        saveAs: { "LiftOver-${sample_id}.vcf.gz" }

    input:
        tuple val(sample_id),
            path(vcf, stageAs: 'inputs/*'),
            path(index, stageAs: 'inputs/*')
        path (header_contigs)
        path (chain_file)

    output:
        tuple val(sample_id), path('liftover.vcf.gz'), emit: liftover_vcf

    script:
        """
        Rscript "${moduleDir}/scripts/liftover-Delly2-vcf.R" \
            --input-vcf "${vcf}" \
            --header-contigs "${header_contigs}" \
            --chain-file "${chain_file}" \
            --output "liftover.vcf.gz"
        """

    stub:
    """
    touch "liftover.vcf.gz"
    """
}

process run_sort_BCFtools {
    container params.docker_image_bcftools

    publishDir path: "${params.output_dir_base}/intermediate/${task.process.replace(':', '/')}",
        pattern: "sorted.vcf.gz",
        mode: "copy",
        enabled: params.save_intermediate_files,
        saveAs: { "BCFtools-sorted-${sample_id}.vcf.gz" }

    input:
        tuple val(sample_id), path(vcf, stageAs: 'inputs/*')

    output:
        tuple val(sample_id), path('sorted.vcf.gz'), emit: sorted_vcf

    script:
        """
        bcftools sort \
            --output-type z \
            --output sorted.vcf.gz \
            ${vcf}
        """

    stub:
        """
        touch sorted.vcf.gz
        """
}

process run_intersect_gnomad {
    container params.docker_image_stablelift

    publishDir path: "${params.output_dir_base}/intermediate/${task.process.replace(':', '/')}",
        pattern: "annotations.Rds",
        mode: "copy",
        enabled: params.save_intermediate_files,
        saveAs: { "LiftOver-${sample_id}-${variant_caller}.Rds" }

    input:
        tuple val(sample_id), path(vcf, stageAs: 'inputs/*')
        path (gnomad_rds)
        val (variant_caller)

    output:
        tuple val(sample_id), path('annotations.Rds'), emit: r_annotations

    script:
        """
        Rscript "${moduleDir}/scripts/extract-vcf-features-SV.R" \
            --variant-caller "${variant_caller}" \
            --input-vcf "${vcf}" \
            --output-rds "annotations.Rds" \
            --gnomad-rds ${gnomad_rds}
        """

    stub:
    """
    touch "annotations.Rds"
    """
}

workflow workflow_extract_sv_annotations {
    take:
    vcf_with_sample_id
    header_contigs
    gnomad_rds
    chain_file
    variant_caller

    main:

    // Step 1: Liftover
    liftover_SV_StableLift(
        vcf_with_sample_id,
        header_contigs,
        chain_file
    )
    run_sort_BCFtools(
        liftover_SV_StableLift.out.liftover_vcf
    )

    // Step 2: Extract features
    run_intersect_gnomad(
        run_sort_BCFtools.out.sorted_vcf,
        gnomad_rds,
        variant_caller
    )

    emit:
    liftover_vcf = run_sort_BCFtools.out.sorted_vcf
    r_annotations = run_intersect_gnomad.out.r_annotations
}
