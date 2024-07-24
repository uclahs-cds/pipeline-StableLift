
process run_sv_liftover{
    container params.docker_image_stablelift

    publishDir path: "${intermediate_filepath}",
        pattern: "liftover.vcf.gz",
        mode: "copy",
        enabled: params.save_intermediate_files,
        saveAs: { "${slug}.vcf.gz" }

    input:
        tuple val(sample_id),
            path(vcf, stageAs: 'inputs/*'),
            path(index, stageAs: 'inputs/*')
        path (header_contigs)
        path (chain_file)

    output:
        tuple val(sample_id), path('liftover.vcf.gz'), emit: liftover_vcf

    script:
        // FIXME Use a more standard path
        intermediate_path = "${params.output_dir_base}/StableLift-${params.stable_version}/intermediate/${task.process}"

        slug = "LiftOver-${sample_id}"
        
        """
        Rscript "${moduleDir}/scripts/liftover-Delly2-vcf.R \
            --input-vcf ${vcf} \
            --header-contigs ${header_contigs} \
            --chain-file ${chain_file} \
            --output "liftover.vcf.gz"
        """
}

process run_intersect_gnomad {
    container params.docker_image_stablelift

    publishDir path: "${intermediate_filepath}",
        pattern: "annotations.Rds",
        mode: "copy",
        enabled: params.save_intermediate_files,
        saveAs: { "${slug}.Rds" }

    input:
        tuple val(sample_id), path(vcf, stageAs: 'inputs/*')
        path (gnomad_rds)
        val (variant_caller)

    output:
        tuple val(sample_id), path('annotations.Rds'), emit: r_annotations

    script:
        // FIXME Use a more standard path
        intermediate_path = "${params.output_dir_base}/StableLift-${params.stable_version}/intermediate/${task.process}"

        slug = "LiftOver-${sample_id}-${variant_caller}"
        
        """
        Rscript ${moduleDir}/scripts/publish/extract-vcf-features-SV.R \
            --variant-caller "${variant_caller}" \
            --input-vcf "${vcf}" \
            --output-rds "annotations.Rds" \
            --gnomad-rds ${gnomad_rds}
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

    run_sv_liftover(
        vcf_with_sample_id,
        header_contigs,
        chain_file
    )

    run_intersect_gnomad(
        run_sv_liftover.out.liftover_vcf,
        gnomad_rds,
        variant_caller
    )

    emit:
    r_annotations = run_intersect_gnomad.out.r_annotations
}
