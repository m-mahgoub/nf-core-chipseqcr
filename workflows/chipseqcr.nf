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

// Check mandatory parameters

// Check Input Sample Sheet
if (params.input) { ch_input = file(params.input) }
else { exit 1, 'Input samplesheet not specified!' }

// Check either Bowtie2 index or Genome fasta files are provided
if (params.bowtie2) { bowtie2index_user =  params.bowtie2}
else if (params.fasta) { ch_fasta = file(params.fasta) }
else { exit 1, "Neither Bowtie2 Index nor Genome Fasta are specified! Provide at least one of these options via '--fasta genome.fa' or '--bowtie2 path/to/bowtie2/index' or via a detectable config file." }

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
include { INPUT_CHECK } from '../subworkflows/local/input_check'

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
include { SAMTOOLS_INDEX               } from '../modules/nf-core/modules/samtools/index/main'

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

    // If Bowtie2 Index is not provides, Build it from fasta:
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
    // Channel Operation: Make A channel with each sample's bam/bai matched to its control's bam/bai
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
    .set { ch_sample_control_bam_bai }
    // Credit for the above two channel operations for nf-core chipseq
    // https://github.com/nf-core/chipseq/blob/813cb384f51901202ff1e937846a6ece4dd199b3/workflows/chipseq.nf

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

    MULTIQC (
        ch_multiqc_files.collect()
    )
    multiqc_report = MULTIQC.out.report.toList()
    ch_versions    = ch_versions.mix(MULTIQC.out.versions)

    // Emit for testing purpose
    emit:ch_sample_control_bam_bai
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
