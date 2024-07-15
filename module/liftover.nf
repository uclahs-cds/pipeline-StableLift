/*
*   Module/process description here
*
*   @input  <name>  <type>  <description>
*   @params <name>  <type>  <description>
*   @output <name>  <type>  <description>
*/

// include { generate_standard_filename } from '../external/pipeline-Nextflow-module/modules/common/generate_standardized_filename/main.nf'

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
