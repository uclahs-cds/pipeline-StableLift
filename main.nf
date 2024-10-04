#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Include processes and workflows here
include { run_validate_PipeVal_with_metadata } from './external/pipeline-Nextflow-module/modules/PipeVal/validate/main.nf' addParams(
    options: [
        docker_image_version: params.pipeval_version,
        main_process: "./" //Save logs in <log_dir>/process-log/run_validate_PipeVal
    ]
)

include { workflow_extract_sv_annotations } from './module/sv_workflow.nf'
include { workflow_extract_snv_annotations } from './module/snv_workflow.nf'
include { workflow_predict_stability } from './module/predict_stability.nf'

// Log info here
log.info """\
    =====================================
    S T A B L E L I F T   P I P E L I N E
    =====================================
    Boutros Lab

    Current Configuration:
    - pipeline:
        name: ${workflow.manifest.name}
        version: ${workflow.manifest.version}

    - input:
        dataset_id: ${params.dataset_id}

        liftover_direction: ${params.liftover_direction}

        variant_caller: ${params.variant_caller}
        rf_model: ${params.rf_model}

        src_fasta_ref:   ${params.src_fasta_ref}
        src_fasta_fai:   ${params.src_fasta_fai}
        src_fasta_dict:  ${params.src_fasta_dict}

        dest_fasta_ref:  ${params.dest_fasta_ref}
        dest_fasta_fai:  ${params.dest_fasta_fai}
        dest_fasta_dict: ${params.dest_fasta_dict}

        chain_file: ${params.chain_file}

    - SV only:
        header_contigs: ${params.getOrDefault('header_contigs', null)}
        gnomad_rds: ${params.getOrDefault('gnomad_rds', null)}

    - SNV only:
        repeat_bed: ${params.getOrDefault('repeat_bed', null)}
        funcotator_data_source: ${params.getOrDefault('funcotator_data_source', null)}

    - output:
        output_dir_base: ${params.output_dir_base}

    - options:
        blcds_registered_dataset: ${params.blcds_registered_dataset}
        ucla_cds: ${params.ucla_cds}

        min_cpus: ${params.min_cpus}
        max_cpus: ${params.max_cpus}

        min_memory: ${params.min_memory}
        max_memory: ${params.max_memory}

    Tools Used:
        BCFtools:   ${params.docker_image_bcftools}
        BEDtools:   ${params.docker_image_bedtools}
        PipeVal:    ${params.docker_image_pipeval}
        SAMTools:   ${params.docker_image_samtools}
        StableLift: ${params.docker_image_stablelift}
        GATK:       ${params.docker_image_gatk}

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
        params.src_fasta_id,
        params.src_fasta_ref,
        params.src_fasta_fai,
        params.src_fasta_dict,
    ] )
    .set { input_ch_src_sequence }

Channel
    .value( [
        params.dest_fasta_id,
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
        .tap { validated_vcf_with_index }
        .map { [it.sample_id, it.vcf, it.index] }
        .set { validated_vcf_tuple }

    // The values of validated_vcf_with_index are maps with keys vcf, index, and sample_id.
    // The values of validated_vcf_tuple are tuples of (sample_id, vcf, index).

    if (params.variant_caller == "Delly2-gSV" || params.variant_caller == "Delly2-sSV") {
        // Take the SV branch
        workflow_extract_sv_annotations(
            validated_vcf_tuple,
            input_ch_src_sequence,
            Channel.value(params.header_contigs),
            Channel.value(params.gnomad_rds),
            Channel.value(params.chain_file)
        )

        workflow_extract_sv_annotations.out.liftover_vcf.set { liftover_vcf }
        workflow_extract_sv_annotations.out.r_annotations.set { r_annotations }

    } else {
        // Take the SNV branch
        workflow_extract_snv_annotations(
            validated_vcf_tuple,
            input_ch_src_sequence,
            input_ch_dest_sequence,
            Channel.value(params.chain_file),
            Channel.value(params.variant_caller)
        )

        workflow_extract_snv_annotations.out.liftover_vcf.set { liftover_vcf }
        workflow_extract_snv_annotations.out.r_annotations.set { r_annotations }
    }

    // Predict stability and apply annotate lifted VCF
    workflow_predict_stability(
        liftover_vcf,
        r_annotations,
        Channel.value(params.rf_model),
        Channel.value(params.variant_caller)
    )
}
