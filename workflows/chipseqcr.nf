/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowChipseqcr.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check YAML file for Heatmaps Blueprint if params.custom_heatmap = true
if (params.custom_heatmap) { file(params.heatmap_blueprint, checkIfExists: true) }

// Check mandatory parameters

// Check Input Sample Sheet
if (params.input) { ch_input = file(params.input) }
else { exit 1, 'Input samplesheet not specified!' }

// Check either Bowtie2 index or Genome fasta files are provided
if (params.bowtie2) { bowtie2index_user =  params.bowtie2}
else if (params.fasta) { ch_fasta = file(params.fasta) }
else { exit 1, "Neither Bowtie2 Index nor Genome Fasta are specified! Provide at least one of these options via '--fasta genome.fa' or '--bowtie2 path/to/bowtie2/index' or via a detectable config file." }

// Construct Channels from Config files
if (params.custom_heatmap) { ch_blueprint = file(params.heatmap_blueprint) }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK             } from '../subworkflows/local/input_check'
include { HEATMAP_BLUEPRINT       } from '../modules/local/heatmap_blueprint'
include { DEEPTOOLS_COMPUTEMATRIX } from '../modules/local/deeptools_computematrix'
include { DEEPTOOLS_PLOTHEATMAP   } from '../modules/local/deeptools_plotheatmap'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                      } from '../modules/nf-core/modules/fastqc/main'
include { MULTIQC                     } from '../modules/nf-core/modules/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'
include { BOWTIE2_ALIGN               } from '../modules/nf-core/modules/bowtie2/align/main'
include { BOWTIE2_BUILD               } from '../modules/nf-core/modules/bowtie2/build/main'
include { SAMTOOLS_SORT               } from '../modules/nf-core/modules/samtools/sort/main'
include { SAMTOOLS_INDEX              } from '../modules/nf-core/modules/samtools/index/main'
include { MACS2_CALLPEAK              } from '../modules/nf-core/modules/macs2/callpeak/main'
include { DEEPTOOLS_BAMCOVERAGE       } from '../modules/nf-core/modules/deeptools/bamcoverage/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow CHIPSEQCR {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        INPUT_CHECK.out.reads,
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // MODULE: Run Bowtie2 Alignment
    //

    // If Bowtie2 Index is not provided, build it from fasta:
    if (!params.bowtie2 && params.fasta) {
       BOWTIE2_BUILD(ch_fasta)
       bowtie2index_fasta = BOWTIE2_BUILD.out.index
       ch_versions = ch_versions.mix(BOWTIE2_BUILD.out.versions)
    }

    // Set Bowtie2 Index Chnannel based on provided input
    if (params.bowtie2) { ch_index = bowtie2index_user}
    else if (!params.bowtie2 && params.fasta)  { ch_index =  bowtie2index_fasta }

    // Perform Bowtie2 Alignment
    BOWTIE2_ALIGN (
        INPUT_CHECK.out.reads,
        ch_index,
        false  // whether to dump unmapped reads
    )
    ch_versions = ch_versions.mix(BOWTIE2_ALIGN.out.versions.first())

    //
    // MODULE: Run Samtools Sort
    //
    SAMTOOLS_SORT (
        BOWTIE2_ALIGN.out.bam
    )

    //
    // MODULE: Run Samtools Index
    //
    SAMTOOLS_INDEX (
        SAMTOOLS_SORT.out.bam
    )


    //
    // Channel Operation: Make A channel emitting bam and bai index of each sample
    //
    SAMTOOLS_SORT
    .out
    .bam
    .join (SAMTOOLS_INDEX.out.bai, by: [0])
    .set { ch_bam_bai }

    //
    // Channel Operation: Make A channel for control samples only
    //
    SAMTOOLS_SORT
    .out
    .bam
    .join (SAMTOOLS_INDEX.out.bai, by: [0])
    .map {
        meta, bam, bai ->
            meta.control ? null : [ meta.id, [ bam ] , [ bai ] ]
    }
    .set { ch_control_bam_bai }

    //
    // Channel Operation: Make A channel with each IP's bam matched to its control's bam (No bai)
    //
    SAMTOOLS_SORT
    .out
    .bam
    .join (SAMTOOLS_INDEX.out.bai, by: [0])
    .map {
        meta, bam, bai ->
            meta.control ? [ meta.control, meta,  [ bam ] ,  [ bai ] ] : null
    }
    .combine(ch_control_bam_bai, by: 0)
    .map { it -> [ it[1] , it[2] , it[4]] }
    .set { ch_ip_control_bam }

    //
    // Channel Operation: Make A channel with each IP's bam/bai matched to its control's bam/bai
    //
    SAMTOOLS_SORT
    .out
    .bam
    .join (SAMTOOLS_INDEX.out.bai, by: [0])
    .map {
        meta, bam, bai ->
            meta.control ? [ meta.control, meta,  [ bam ] ,  [ bai ] ] : null
    }
    .combine(ch_control_bam_bai, by: 0)
    .map { it -> [ it[1] , it[2] + it[4], it[3] + it[5] ] }
    .set { ch_ip_control_bam_bai }
    // Credit for the code of the above two channel operations for nf-core chipseq
    // https://github.com/nf-core/chipseq/blob/813cb384f51901202ff1e937846a6ece4dd199b3/workflows/chipseq.nf

    //
    // MODULE: Run MACS2 Peak Calling
    //
    MACS2_CALLPEAK (
        ch_ip_control_bam,
        params.gsize
    )
    ch_versions = ch_versions.mix(MACS2_CALLPEAK.out.versions.first())

    //
    // MODULE: Run Deeptools bamCoverage
    //
    DEEPTOOLS_BAMCOVERAGE (
        ch_bam_bai
    )
    ch_versions = ch_versions.mix(DEEPTOOLS_BAMCOVERAGE.out.versions.first())

    //
    // Optional Subworkflow: if custom Heatmap is required
    //
    if (params.custom_heatmap) {
        HEATMAP_BLUEPRINT (
            ch_blueprint
        )
        ch_versions = ch_versions.mix(HEATMAP_BLUEPRINT.out.versions)

        //
        // Channel Operation: Make A channel that collects all bed files for MACS2_CALLPEAK in one channel
        //
        MACS2_CALLPEAK
        .out
        .peak
        .map { it -> [it[1]]}
        .collect()
        .set { ch_bed_all }

        //
        // Channel Operation: Make A channel that collects all bigwig files for DEEPTOOLS_BAMCOVERAGE in one channel
        //
        DEEPTOOLS_BAMCOVERAGE
        .out
        .bigwig
        .map { it -> [it[1]]}
        .collect()
        .set { ch_bigwig_all }

        //
        // Channel Operation: Make A channel that emits string of all bed and bigwig files' names
        //                    This is a mix of files staged by path(bed) and path(local_remote_files)
        //
        HEATMAP_BLUEPRINT.out.str_for_deeptool
        .splitCsv(header:true, sep:'\t')
        .map { row -> [ row.plot_name, row.bed_files_inLine, row.bigwig_files_inLine ] }
        .set { ch_files_str_for_deeptool }


        //
        // Channel Operation: Make A channel that emits paths for local or remote files (not generated in the pipline) to be staged
        //
        HEATMAP_BLUEPRINT.out.path_for_deeptool
        .splitCsv(header:false, sep:'\t')
        .collect()
        .set { ch_path_for_deeptool }


        DEEPTOOLS_COMPUTEMATRIX (
            ch_bed_all, //
            ch_bigwig_all, //
            ch_path_for_deeptool, //
            ch_files_str_for_deeptool
        )
        ch_versions = ch_versions.mix(DEEPTOOLS_COMPUTEMATRIX.out.versions)

        //
        // Channel Operation: Make A channel that emits string for labels of regions and samples, matched to matrix files
        // 
        
        // (1) channel that emits string for labels of regions and samples only
        HEATMAP_BLUEPRINT.out.str_for_deeptool
        .splitCsv(header:true, sep:'\t')
        .map { row -> [ row.plot_name, row.bed_labels_inLine_with_quotes, row.bigwig_labels_inLine_with_quotes ] }
        .set { ch_labels_for_deeptool }
        // (2) channel that emits plot name and assocaited matrix file
        DEEPTOOLS_COMPUTEMATRIX.out.matrix
        .map { it -> [ it.getSimpleName(), it ] }
        .set { ch_matrix }
        // (3) Join the two channels
        ch_matrix
        .join( ch_labels_for_deeptool )
        .set { ch_heatmap_input }

        DEEPTOOLS_PLOTHEATMAP (
            ch_heatmap_input
        )
        ch_versions = ch_versions.mix(DEEPTOOLS_PLOTHEATMAP.out.versions)
    }

    // #####################################################################################################################################
    // #####################################################################################################################################
    // #####################################################################################################################################
    // #####################################################################################################################################

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowChipseqcr.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(BOWTIE2_ALIGN.out.log.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(MACS2_CALLPEAK.out.xls.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect()
    )
    multiqc_report = MULTIQC.out.report.toList()
    ch_versions    = ch_versions.mix(MULTIQC.out.versions)

    // #####################################################################################################################################
    // #####################################################################################################################################
    // #####################################################################################################################################
    // #####################################################################################################################################



    // Emit for testing purpose
    // emit: ch_bed_all
    // emit: ch_bam_bai
    // emit: ch_labels_for_deeptool
    // emit: ch_files_str_for_deeptool
    // emit: ch_heatmap_input
    emit: Channel.empty()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
