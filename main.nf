#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Include processes and workflows here
include { run_validate_PipeVal_with_metadata } from './external/pipeline-Nextflow-module/modules/PipeVal/validate/main.nf' addParams(
    options: [
        docker_image_version: params.pipeval_version,
        main_process: "./" //Save logs in <log_dir>/process-log/run_validate_PipeVal
    ]
)

include { raw_liftover } from './module/liftover.nf'
include { run_funcotator } from './module/funcotator.nf'
include { apply_annotations } from './module/annotations.nf'
include { extract_features} from './module/extract_features.nf'

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

Channel
    .value( [
        params.funcotator_data.src_reference_id,
        params.src_fasta_ref,
        params.src_fasta_fai,
        params.src_fasta_dict,
    ] )
    .set { input_ch_src_sequence }

Channel
    .value( [
        params.funcotator_data.dest_reference_id,
        params.dest_fasta_ref,
        params.dest_fasta_fai,
        params.dest_fasta_dict
    ] )
    .set { input_ch_dest_sequence }

// Main workflow here
workflow {

    // Currently this is written for a single sample_id and VCF file, but
    // abstract that away
    Channel.of ([
            vcf: params.input.vcf,
            index: indexFile(params.input.vcf),
            sample_id: params.sample_id
        ]).set { vcf_with_index }

    // The values of vcf_with_index are maps with keys vcf, index, and sample_id.

    // Run the input VCF and TBI files through PipeVal
    vcf_with_index
        .flatMap { sample ->
            [
                [sample.vcf, [[sample_id: sample.sample_id], "vcf"]],
                [sample.index, [[sample_id: sample.sample_id], "index"]]
            ]
        } | run_validate_PipeVal_with_metadata

    // Save the validation result
    run_validate_PipeVal_with_metadata.out.validation_result
        .collectFile(
            name: 'input_validation.txt',
            newLine: true,
            storeDir: "${params.output_dir_base}/validation"
        )

    run_validate_PipeVal_with_metadata.out.validated_file
        .map { filename, metadata -> [metadata[0].sample_id, metadata[0] + [(metadata[1]): filename]] }
        .groupTuple()
        .map { it[1].inject([:]) { result, i -> result + i } }
        .set { validated_vcf_with_index }

    // The values of validated_vcf_with_index are maps with keys vcf, index, and sample_id.
    raw_liftover(
        validated_vcf_with_index.map { [it.sample_id, it.vcf, it.index] },
        input_ch_src_sequence,
        input_ch_dest_sequence,
        Channel.value(params.chain_file)
    )

    run_funcotator(
        raw_liftover.out.liftover_vcf_with_index,
        input_ch_dest_sequence,
        Channel.value(params.funcotator_data.data_source)
    )

    apply_annotations(
        run_funcotator.out.funcotator_vcf,
        input_ch_dest_sequence
    )

    extract_features(
        apply_annotations.out.annotated_vcf
    )
}
