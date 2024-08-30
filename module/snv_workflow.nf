include { workflow_apply_snv_annotations } from './snv_annotations.nf'

process run_liftover_BCFtools {
    container params.docker_image_bcftools

    publishDir path: "${params.output_dir_base}/intermediate/${task.process.replace(':', '/')}",
        pattern: "{reject,liftover}.vcf.gz{,.tbi}",
        mode: "copy",
        enabled: params.save_intermediate_files,
        saveAs: { filename -> "LiftOver-${sample_id}-${src_fasta_id}-to-${dest_fasta_id}-${filename}" }

    input:
        tuple val(sample_id), path(vcf), path(index)
        tuple val(src_fasta_id), path(src_fasta_ref), path(src_fasta_fai), path(src_fasta_dict)
        tuple val(dest_fasta_id), path(dest_fasta_ref), path(dest_fasta_fai), path(dest_fasta_dict)
        path (chain_file)

    output:
        tuple val(sample_id), path('liftover.vcf.gz'), path('liftover.vcf.gz.tbi'), emit: liftover_vcf_with_index

    script:
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

    stub:
    """
    touch "liftover.vcf.gz"
    touch "liftover.vcf.gz.tbi"
    touch "reject.vcf.gz"
    """
}

process extract_VCF_features_StableLift {
    container params.docker_image_stablelift
    containerOptions "-v ${moduleDir}:${moduleDir}"

    publishDir path: "${params.output_dir_base}/intermediate/${task.process.replace(':', '/')}",
        pattern: "features.Rds",
        mode: "copy",
        enabled: params.save_intermediate_files,
        saveAs: { "StableLift-${sample_id}.Rds" }

    input:
    tuple val(sample_id), path(vcf), path(index)

    output:
    tuple val(sample_id), path('features.Rds'), emit: r_annotations

    script:
    """
    Rscript "${moduleDir}/scripts/extract-VCF-features.R" \
        --input-vcf "${vcf}" \
        --variant-caller ${params.variant_caller} \
        --output-rds "features.Rds"
    """

    stub:
    """
    touch "features.Rds"
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

    // We want to do all of the annotating with the GRCh38 / hg38 reference. If
    // the liftover is going from h38 to hg19, defer until after annotations
    if (params.liftover_direction == "GRCh37ToGRCh38") {
        // Step 1: Liftover
        run_liftover_BCFtools(
            vcf_with_sample_id,
            src_sequence,
            dest_sequence,
            chain_file
        )

        // Step 2: Annotate with GRCh38
        workflow_apply_snv_annotations(
            run_liftover_BCFtools.out.liftover_vcf_with_index,
            dest_sequence
        )

        workflow_apply_snv_annotations.out.annotated_vcf.set { annotated_vcf_with_index }

    } else {
        // Step 1: Annotate with GRCh38
        workflow_apply_snv_annotations(
            vcf_with_sample_id,
            src_sequence
        )

        // Step 2: Liftover
        run_liftover_BCFtools(
            workflow_apply_snv_annotations.out.annotated_vcf,
            src_sequence,
            dest_sequence,
            chain_file
        )

        run_liftover_BCFtools.out.liftover_vcf_with_index.set { annotated_vcf_with_index }
    }

    // Step 3: Extract features
    // FIXME Parallelize HaplotypeCaller
    extract_VCF_features_StableLift(
        annotated_vcf_with_index
    )

    // For consistency with the SV branch, remove the index file from the
    // output VCF channel
    annotated_vcf_with_index
        .map { sample_id, vcf, index -> [sample_id, vcf] }
        .set { annotated_vcf }

    emit:
    liftover_vcf = annotated_vcf
    r_annotations = extract_VCF_features_StableLift.out.r_annotations
}

