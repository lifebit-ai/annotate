#!/usr/bin/env nextflow
/*
========================================================================================
                         lifebit-ai/annotate
========================================================================================
 lifebit-ai/annotate Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/lifebit-ai/annotate
----------------------------------------------------------------------------------------
*/

def helpMessage() {
    // TODO nf-core: Add to this help message with new command line parameters
    log.info nfcoreHeader()
    log.info"""

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run main.mf ...

    Mandatory arguments:
                                    Available: conda, docker, singularity, test, awsbatch, <institute> and more

    Options:

    Other options:
      --outdir [file]                 The output directory where the results will be saved
      -name [str]                     Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic
    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}


// TODO nf-core: Add any reference files that are needed

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
    custom_runName = workflow.runName
}

ch_output_docs = file("$baseDir/docs/output.md", checkIfExists: true)


// Header log info
log.info nfcoreHeader()
def summary = [:]
if (workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name']         = custom_runName ?: workflow.runName
// TODO nf-core: Report custom parameters here
summary['Max Resources']    = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
if (workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['Output dir']       = params.outdir
summary['Launch dir']       = workflow.launchDir
summary['Working dir']      = workflow.workDir
summary['Script dir']       = workflow.projectDir
summary['User']             = workflow.userName
summary['Config Profile'] = workflow.profile
if (params.config_profile_description) summary['Config Description'] = params.config_profile_description
if (params.config_profile_contact)     summary['Config Contact']     = params.config_profile_contact
if (params.config_profile_url)         summary['Config URL']         = params.config_profile_url
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "-\033[2m--------------------------------------------------\033[0m-"


Channel.from(summary.collect{ [it.key, it.value] })
    .map { k,v -> "<dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }
    .reduce { a, b -> return [a, b].join("\n            ") }
    .map { x -> """
    id: 'nf-core-annotate-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'lifebit-ai/annotate Workflow Summary'
    section_href: 'https://github.com/lifebit-ai/annotate'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
            $x
        </dl>
    """.stripIndent() }
    .set { ch_workflow_summary }



// Defining input channels

Channel.fromPath(params.predicted_ancestries)
    .set { ch_infer_ancestry }

Channel.fromPath(params.unrelated_list)
    .set { ch_unrelated_list }



/*
 * STEP 1 - prep_hwe
 */

process prep_hwe {
    publishDir "${params.outdir}/prep_hwe_pop_files", mode: params.publish_dir_mode

    input:
    file(predicted_ancestries) from ch_infer_ancestry
    file(unrelated_list) from ch_unrelated_list

    output:
    file("*_unrelated_pop.keep") into ch_unrelated_by_pop_keep_files

    script:
    """
    prep_hwe.R --predicted_ancestries='${predicted_ancestries}' \
               --unrelated_list='${unrelated_list}'
    """
}

ch_unrelated_by_pop_keep_files.flatten().view()

/*
 * Completion notification
 */
workflow.onComplete {

    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";

    if (workflow.stats.ignoredCount > 0 && workflow.success) {
        log.info "-${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}-"
        log.info "-${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCount} ${c_reset}-"
        log.info "-${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCount} ${c_reset}-"
    }

    if (workflow.success) {
        log.info "-${c_purple}[lifebit-ai/annotate]${c_green} Pipeline completed successfully${c_reset}-"
    } else {
        log.info "-${c_purple}[lifebit-ai/annotate]${c_red} Pipeline completed with errors${c_reset}-"
    }

}


def nfcoreHeader() {
    // Log colors ANSI codes
    c_black = params.monochrome_logs ? '' : "\033[0;30m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_white = params.monochrome_logs ? '' : "\033[0;37m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";

    return """    -${c_dim}--------------------------------------------------${c_reset}-
                                            ${c_green},--.${c_black}/${c_green},-.${c_reset}
    ${c_blue}        ___     __   __   __   ___     ${c_green}/,-._.--~\'${c_reset}
    ${c_blue}  |\\ | |__  __ /  ` /  \\ |__) |__         ${c_yellow}}  {${c_reset}
    ${c_blue}  | \\| |       \\__, \\__/ |  \\ |___     ${c_green}\\`-._,-`-,${c_reset}
                                            ${c_green}`._,._,\'${c_reset}
    ${c_purple}  lifebit-ai/annotate v${workflow.manifest.version}${c_reset}
    -${c_dim}--------------------------------------------------${c_reset}-
    """.stripIndent()
}
