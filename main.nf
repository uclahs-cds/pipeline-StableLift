#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Include processes and workflows here
include { run_validate_PipeVal } from './external/pipeline-Nextflow-module/modules/PipeVal/validate/main.nf' addParams(
    options: [
        docker_image_version: params.pipeval_version,
        main_process: "./" //Save logs in <log_dir>/process-log/run_validate_PipeVal
    ]
)
include { generate_standard_filename } from './external/pipeline-Nextflow-module/modules/common/generate_standardized_filename/main.nf'
include { tool_name_command_name } from './module/module-name' addParams(
    workflow_output_dir: "${params.output_dir_base}/ToolName-${params.toolname_version}",
    workflow_log_output_dir: "${params.log_output_dir}/process-log/ToolName-${params.toolname_version}"
    )

// Log info here
log.info """\
        ======================================
        T E M P L A T E - N F  P I P E L I N E
        ======================================
        Boutros Lab

        Current Configuration:
        - pipeline:
            name: ${workflow.manifest.name}
            version: ${workflow.manifest.version}

        - input:
            input a: ${params.variable_name}
            ...

        - output:
            output a: ${params.output_path}
            ...

        - options:
            option a: ${params.option_name}
            ...

        Tools Used:
            tool a: ${params.docker_image_name}

        ------------------------------------
        Starting workflow...
        ------------------------------------
        """
        .stripIndent()

def indexFile(bam_or_vcf) {
    if(bam_or_vcf.endsWith('.bam')) {
        return "${bam_or_vcf}.bai"
        }
    else if(bam_or_vcf.endsWith('vcf.gz')) {
        return "${bam_or_vcf}.tbi"
        }
    else {
        throw new Exception("Index file for ${bam_or_vcf} file type not supported. Use .bam or .vcf.gz files.")
        }
    }

// Channels here
Channel
    .fromList(params.samples_to_process)
    .map { sample ->
        return tuple(sample.orig_id, sample.id, sample.path, sample.sample_type)
    }
    .set { samplesToProcessChannel }

Channel
    .fromList(params.samples_to_process)
    .map{ it -> [it['path'], indexFile(it['path'])] }
    .flatten()
    .set { files_to_validate_ch }

// Possible reference channel
Channel
    .from(
        params.reference,
        params.reference_index,
        params.reference_dict
        )
    .set { reference_ch }

// Decription of input channel
Channel
    .fromPath(params.variable_name)
    .ifEmpty { error "Cannot find: ${params.variable_name}" }
    .set { input_ch_variable_name }

files_to_validate_ch = files_to_validate_ch
    .mix(reference_ch)
    .mix(input_ch_variable_name)

// Main workflow here
workflow {
    // Validation process
    run_validate_PipeVal(
        files_to_validate_ch
        )

    run_validate_PipeVal.out.val_file.collectFile(
        name: 'input_validation.txt', newLine: true,
        storeDir: "${params.output_dir_base}/validation"
        )

    // Workflow or process
    tool_name_command_name(
        samplesToProcessChannel,
        input_ch_variable_name
        )
}
