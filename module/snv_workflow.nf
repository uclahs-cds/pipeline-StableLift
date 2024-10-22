include { workflow_annotate_snvs } from './snv_annotations.nf'

process run_liftover_BCFtools {
    container params.docker_image_bcftools

    publishDir path: "${params.output_dir_base}/intermediate/${task.process.replace(':', '/')}",
        pattern: "{reject,liftover}.vcf.gz{,.tbi}",
        mode: "copy",
        enabled: params.save_intermediate_files,
        saveAs: { filename -> "${sample_id}_LiftOver-${dest_fasta_id}_${filename}" }

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
            --af-tags "" \
            --reject-type z \
            --reject "reject.vcf.gz" | \
            bcftools view \
                --output-type u \
                -e 'REF="." | ALT="." | POS<10' | \
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

process split_vcf_for_feature_extraction {
    container params.docker_image_bcftools
    
    input:
    tuple val(sample_id), path(vcf), path(index)

    output:
    tuple val(sample_id), path('split_vcf'), emit: split_vcf

    script:
    """
    mkdir -p split_vcf
    vcf=\$(readlink -f "${vcf}")
    if [ \$(du -m "\${vcf}" | cut -f1) -gt 1024 ]; then
        bcftools view -H "\${vcf}" |
            split -C 1G -d -a2 \
                --additional-suffix=".vcf" \
                - "split_vcf/${sample_id}_split_"
        for file in split_vcf/*; do
            bcftools view -h "\${vcf}" > "\${file}~"
            cat "\${file}" >> "\${file}~"
            mv "\${file}~" "\${file}"
        done
    else
        ln -s \${vcf} "split_vcf/${sample_id}.vcf"
    fi
    """
}

process extract_VCF_features_StableLift {
    container params.docker_image_stablelift
    containerOptions "-v ${moduleDir}:${moduleDir}"

    cpus { params.getOrDefault('extract_features_cpus', 4) }

    publishDir path: "${params.output_dir_base}/intermediate/${task.process.replace(':', '/')}",
        pattern: "features.Rds",
        mode: "copy",
        enabled: params.save_intermediate_files,
        saveAs: { "StableLift-${sample_id}.Rds" }

    input:
    tuple val(sample_id), path(vcf_dir)
    val variant_caller

    output:
    tuple val(sample_id), path('features.Rds'), emit: r_annotations

    script:
    """
    Rscript "${moduleDir}/scripts/extract-VCF-features.R" \
        --input-dir "${vcf_dir}" \
        --output-rds "features.Rds" \
        --ncore ${task.cpus} \
        --variant-caller ${variant_caller}
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

    // Step 1: LiftOver
    run_liftover_BCFtools(
        vcf_with_sample_id,
        src_sequence,
        dest_sequence,
        chain_file
    )

    // Step 2: Annotate
    workflow_annotate_snvs(
        run_liftover_BCFtools.out.liftover_vcf_with_index,
        dest_sequence
    )

    // Step 3: Split large VCFs
    split_vcf_for_feature_extraction(
        workflow_annotate_snvs.out.annotated_vcf
    )

    // Step 4: Extract features
    extract_VCF_features_StableLift(
        split_vcf_for_feature_extraction.out.split_vcf,
        variant_caller
    )

    // For consistency with the SV branch, remove the index file from the output VCF channel
    workflow_annotate_snvs.out.annotated_vcf
        .map { sample_id, vcf, index -> [sample_id, vcf] }
        .set { annotated_vcf }

    emit:
    liftover_vcf = annotated_vcf
    r_annotations = extract_VCF_features_StableLift.out.r_annotations
}
