/*
*   Module/process description here
*
*   @input  <name>  <type>  <description>
*   @params <name>  <type>  <description>
*   @output <name>  <type>  <description>
*/

include { generate_standard_filename } from '../external/pipeline-Nextflow-module/modules/common/generate_standardized_filename/main.nf'

process run_Funcotator_GATK {
    container params.docker_image_gatk

    publishDir path: "${intermediate_filepath}",
        pattern: "output.vcf.gz",
        mode: "copy",
        enabled: params.save_intermediate_files,
        saveAs: { "${slug}.vcf.gz" }

    input:
        tuple val(sample_id),
            path(vcf, stageAs: 'inputs/*'),
            path(index, stageAs: 'inputs/*')
        tuple val(dest_fasta_id), path(dest_fasta_ref), path(dest_fasta_fai), path(dest_fasta_dict)
        path (funcotator_sources)

    output:
        tuple val(sample_id), path('output.vcf.gz'), emit: funcotator_vcf

    script:
        intermediate_filepath = "${params.output_dir_base}/GATK-${params.gatk_version}/intermediate/${task.process}"

        slug = "Funcotator-${sample_id}-${dest_fasta_id}"
        
        """
        gatk Funcotator \
            --variant "${vcf}" \
            --reference "${dest_fasta_ref}" \
            --ref-version "${dest_fasta_id}" \
            --data-sources-path "${funcotator_sources}" \
            --output-file-format VCF \
            --output "output.vcf.gz"
        """
}
