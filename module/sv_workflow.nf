process run_SV_liftover_and_annotate {
    container params.docker_image_stablelift

    publishDir path: "${params.output_dir_base}/intermediate/${task.process.replace(':', '/')}",
        pattern: "liftover.{vcf.gz,Rds}",
        mode: "copy",
        enabled: params.save_intermediate_files,
        saveAs: { "LiftOver-${sample_id}.vcf.gz" }

    input:
        tuple val(sample_id), path(vcf, stageAs: 'inputs/*')
        val(source_grch_label)
        path (header_contigs)
        path (chain_file)
        path (gnomad_rds)

    output:
        tuple val(sample_id), path('annotations.vcf.gz'), emit: liftover_vcf
        tuple val(sample_id), path('annotations.Rds'), emit: r_annotations

    script:
        """
        Rscript "${moduleDir}/scripts/extract-VCF-features-SV.R" \
            --input-vcf "${vcf}" \
            --source-build "${source_grch_label}" \
            --chain-file "${chain_file}" \
            --header-contigs "${header_contigs}" \
            --gnomad-rds ${gnomad_rds} \
            --output-vcf "annotations.vcf.gz" \
            --output-rds "annotations.Rds"
        """

    stub:
    """
    touch "annotations.Rds"
    touch "annotations.vcf.gz"
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

workflow workflow_extract_sv_annotations {
    take:
    vcf_with_sample_id
    src_sequence
    header_contigs
    gnomad_rds
    chain_file

    main:

    run_SV_liftover_and_annotate(
        // We don't need the index file
        vcf_with_sample_id.map{ [it[0], it[1]] },

        // We only need the sample ID
        src_sequence.map{ ["hg19": "GRCh37", "hg38": "GRCh38"][it[0]] },

        header_contigs,
        chain_file,
        gnomad_rds
    )

    run_sort_BCFtools(
        run_SV_liftover_and_annotate.out.liftover_vcf
    )

    emit:
    liftover_vcf = run_sort_BCFtools.out.sorted_vcf
    r_annotations = run_SV_liftover_and_annotate.out.r_annotations
}
