include { workflow_apply_snv_annotations } from './module/snv_annotations.nf'

process run_liftover_BCFtools {
    container params.docker_image_bcftools

    publishDir path: "${intermediate_path}",
        pattern: "reject.vcf.gz",
        mode: "copy",
        enabled: params.save_intermediate_files,
        saveAs: { "${slug}-reject.vcf.gz" }

    publishDir path: "${intermediate_path}",
        pattern: "liftover.vcf.gz",
        mode: "copy",
        enabled: params.save_intermediate_files,
        saveAs: { "${slug}.vcf.gz" }

    publishDir path: "${intermediate_path}",
        pattern: "liftover.vcf.gz.tbi",
        mode: "copy",
        enabled: params.save_intermediate_files,
        saveAs: { "${slug}.vcf.gz.tbi" }

    input:
        tuple val(sample_id), path(vcf), path(index)
        tuple val(src_fasta_id), path(src_fasta_ref), path(src_fasta_fai), path(src_fasta_dict)
        tuple val(dest_fasta_id), path(dest_fasta_ref), path(dest_fasta_fai), path(dest_fasta_dict)
        path (chain_file)

    output:
        tuple val(sample_id), path('liftover.vcf.gz'), path('liftover.vcf.gz.tbi'), emit: liftover_vcf_with_index

    script:
        // FIXME Use a more standard path
        intermediate_path = "${params.output_dir_base}/BCFtools-${params.bcftools_version}/intermediate/${task.process}"

        slug = "LiftOver-${sample_id}-${src_fasta_id}-to-${dest_fasta_id}"
        
        """
        bcftools +liftover \
            --output-type u \
            "${vcf}" \
            -- \
            --src-fasta-ref "${src_fasta_ref}" \
            --fasta-ref "${dest_fasta_ref}" \
            --chain "${chain_file}" \
            --reject-type z \
            --reject "reject.vcf.gz" | \
            bcftools view \
                --output-type u \
                -e 'REF="." | ALT="."' | \
            bcftools sort \
                --output-type z \
                --write-index=tbi \
                --output "liftover.vcf.gz"
        """
}

process extract_VCF_features_StableLift {
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

workflow workflow_extract_snv_annotations {
    take:
    vcf_with_sample_id
    src_sequence
    dest_sequence
    chain_file
    variant_caller

    main:

    // Step 1: Liftover
    run_liftover_BCFtools(
        vcf_with_sample_id,
        src_sequence,
        dest_sequence,
        chain_file
    )

    // Step 2: Annotate
    workflow_apply_snv_annotations(
        run_liftover_BCFtools.out.liftover_vcf_with_index,
        dest_sequence
    )

    // Step 3: Extract features
    // FIXME Parallelize HaplotypeCaller
    extract_VCF_features_StableLift(
        workflow_apply_snv_annotations.out.annotated_vcf
    )

    emit:
    liftover_vcf = workflow_apply_snv_annotations.out.annotated_vcf
    r_annotations = extract_VCF_features_StableLift.out.r_annotations
}
