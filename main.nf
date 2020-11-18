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

      To run pipeline with example test data:

        nextflow run main.nf -profile test

      The typical command for running the pipeline is as follows:

        nextflow run main.mf --input input.csv \\
            --predicted_ancestries predicted_ancestries.tsv \\
            --unrelated_list autosomes_LD_pruned_1kgp3Intersect_triangle_HWE1_5.king.cutoff.in.id \\
            --siteqc_results_dir siteqc_results

    Mandatory arguments:
        --input [file]                  File with list of full paths to bcf files and their indexes.
                                        Bcf files can be compressed but in a readable for bcftools format.
                                        Example:
                                        #-----my_bcf_files_list.csv-----------#
                                        | bcf,index                           |
                                        | <file1.bcf>,<file1.bcf.idx>         |
                                        | <file2.bcf.gz>,<file2.bcf.gz.csi>   |
                                        | <file3.bcf.bgz>,<file3.bcf.bgz.tbx> |
                                        #-------------------------------------#
                                        The name of the files must be consistent across files
                                        and follow a specific pattern:
                                        {name}_{CHR}_{START_POS}_{END_POS}.bcf.gz
                                        Example:
                                        test_all_chunks_merged_norm_chr10_53607810_55447336.bcf.gz
                                        Consistency is important here as a variable ('region')
                                        is extracted from the filename.

        --predicted_ancestries [file]   File containing predicted ancestry values for participants.
                                        Expected to be generated by Ancestry and relatedness pipeline,
                                        by process king_coefficients. In that pipeline it is called
                                        predicted_ancestries.tsv

        --unrelated_list [file]         File containing list of unrelated participants (platekeys).
                                        Expected to be generated by Ancestry and relatedness pipeline,
                                        by process king_coefficients. In that pipeline it is called
                                        autosomes_LD_pruned_1kgp3Intersect_triangle_HWE1_5.king.cutoff.in.id

        --siteqc_results_dir [dir]      The results directory of SiteQC pipeline, that has all metrics files
                                        generated from bcf/vcf files provided with --input option.

    Optional arguments:
        -profile [profile]              Avaialbe options: standard (default), test (sufficient minimal test run).

    Other options:
        --outdir [file]                 The output directory where the results will be saved.
        -name [str]                     Name for the pipeline run. If not specified, Nextflow will automatically
                                        generate a random mnemonic.
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
summary['Input bcf list']   = params.input
summary['Predicted ancestries'] = params.predicted_ancestries
summary['Unrelated list']   = params.unrelated_list
summary['SiteQC dir']       = params.siteqc_results_dir
summary['Output dir']       = params.outdir
summary['Launch dir']       = workflow.launchDir
summary['Working dir']      = workflow.workDir
summary['Script dir']       = workflow.projectDir
summary['User']             = workflow.userName
summary['Config Profile']   = workflow.profile
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

/*
 * Check all important required inputs
 */

// Check if user provided input csv file containing paths to bcf files and their indexes
if (!params.input) exit 1, "The list of input bcf/vcf files was not provided. \nPlease specify a csv file containing paths to bcf/vcf files and their indexes with --input [file] option. \nUse --help option for more information."

// Check if user provided input csv file containing paths to bcf files and their indexes
if (!params.predicted_ancestries) exit 1, "The participants' predicted ancestry file was not specified. \nPlease provide it with --predicted_ancestries [file] option. \nUse --help option for more information."

// Check if user provided input csv file containing paths to bcf files and their indexes
if (!params.unrelated_list) exit 1, "The list of unrelated participants was not specified. \nPlease provide it with --unrelated_list [file] option. \nUse --help option for more information."

// Check if user provided SiteQC results directory
if (!params.siteqc_results_dir) exit 1, "The SiteQC results folder was not specified. \nPlease provide it with --siteqc_results_dir [dir] option. \nUse --help option for more information."

// Check if user provided female sample list (xx.txt)
if (!params.xx_sample_ids) exit 1, "The female sample lits was not specified. \nPlease provide it with --xx_sample_ids [file] option. \nUse --help option for more information."


// Shortening the parameters

siteqc_results_dir = params.siteqc_results_dir

startfiles_dir             = params.startfiles_dir
startfiles_prefix          = params.startfiles_prefix
missing1_dir               = params.missing1_dir
missing1_prefix            = params.missing1_prefix
missing2_dir               = params.missing2_dir
missing2_prefix            = params.missing2_prefix
median_coverage_all_dir    = params.median_coverage_all_dir
median_coverage_all_prefix = params.median_coverage_all_prefix
medianNonMiss_dir          = params.medianNonMiss_dir
medianNonMiss_prefix       = params.medianNonMiss_prefix
medianGQ_dir               = params.medianGQ_dir
medianGQ_prefix            = params.medianGQ_prefix
hetAll_dir                 = params.hetAll_dir
hetAll_prefix              = params.hetAll_prefix
hetPass_dir                = params.hetPass_dir
hetPass_prefix             = params.hetPass_prefix
MendErr_dir                = params.MendErr_dir
MendErr_prefix             = params.MendErr_prefix
N_samples_suffix           = params.N_samples_suffix
AC_counts_dir              = params.AC_counts_dir
AC_counts_suffix           = params.AC_counts_suffix


// Defining input channels


Channel.fromPath(params.input)
    .ifEmpty { exit 1, "Input .csv list of input tissues not found at ${params.input}. Is the file path correct?" }
    .splitCsv(sep: ',',  skip: 1)
    .map { bcf, index -> ['chr'+file(bcf).simpleName.split('_chr').last() , file(bcf), file(index)] }
    .into { ch_bcfs_p_hwe; ch_bcfs_make_header; ch_bcfs_final_annotation }


Channel.fromPath(params.predicted_ancestries)
    .set { ch_infer_ancestry }

Channel.fromPath(params.unrelated_list)
    .set { ch_unrelated_list }


// SiteQC results channels:

ch_startfiles = Channel.fromPath("${siteqc_results_dir}/${startfiles_dir}/${startfiles_prefix}*")
    .map { path -> ['chr'+file(path).simpleName.split('chr').last().replaceAll("_X[XY]",""), file(path)] }
    .groupTuple()
    // First item in tuple carries region value - bcf file identity, and is deduced from the file name.
    // For chrX there are special cases where males' and females' X chromosomes have to be treated differently.
    // In such cases for chrX regions there are pairs of metrics files, which have additional _XX _XY suffixes.
    // To deal with such cases, for region value we specifically remove the suffix, and pairs of such _XX and _XY
    // files are grouped in one tuple with .groupTuple() operator by the region value that corresponds to the same bcf.
    // Not all metrics files have such special handlingfor chrX.

ch_miss1 = Channel.fromPath("${siteqc_results_dir}/${missing1_dir}/${missing1_prefix}*")
    .map { path -> ['chr'+file(path).simpleName.split('chr').last().replaceAll("_X[XY]",""), file(path)] }
    .groupTuple()

ch_miss2 = Channel.fromPath("${siteqc_results_dir}/${missing2_dir}/${missing2_prefix}*")
    .map { path -> ['chr'+file(path).simpleName.split('chr').last().replaceAll("_X[XY]",""), file(path)] }
    .groupTuple()

ch_medianCoverageAll = Channel.fromPath("${siteqc_results_dir}/${median_coverage_all_dir}/${median_coverage_all_prefix}*")
    .map { path -> ['chr'+file(path).simpleName.split('chr').last().replaceAll("_X[XY]",""), file(path)] }
    .groupTuple()

ch_medianNonMiss = Channel.fromPath("${siteqc_results_dir}/${medianNonMiss_dir}/${medianNonMiss_prefix}*")
    .map { path -> ['chr'+file(path).simpleName.split('chr').last().replaceAll("_X[XY]",""), file(path)] }
    .groupTuple()

ch_medianGQ = Channel.fromPath("${siteqc_results_dir}/${medianGQ_dir}/${medianGQ_prefix}*")
    .map { path -> ['chr'+file(path).simpleName.split('chr').last().replaceAll("_X[XY]",""), file(path)] }
    .groupTuple()

ch_hetAll = Channel.fromPath("${siteqc_results_dir}/${hetAll_dir}/${hetAll_prefix}*")
    .map { path -> ['chr'+file(path).simpleName.split('chr').last().replaceAll("_XX",""), file(path)] }

ch_hetPass = Channel.fromPath("${siteqc_results_dir}/${hetPass_dir}/${hetPass_prefix}*")
    .map { path -> ['chr'+file(path).simpleName.split('chr').last().replaceAll("_XX",""), file(path)] }

ch_MendErr = Channel.fromPath("${siteqc_results_dir}/${MendErr_dir}/${MendErr_prefix}*")
    .map { path -> ['chr'+file(path).simpleName.split('chr').last(), file(path)] }

ch_AC_counts = Channel.fromPath("${siteqc_results_dir}/${AC_counts_dir}/*${AC_counts_suffix}")
    .map { path -> [file(path).simpleName.split('_AC').first(), file(path)] }

ch_N_samples_autosomes = Channel.fromPath("${siteqc_results_dir}/*${N_samples_suffix}")
    .map { path -> [ "autosmomal_region", file(path)] }
// Must always look as following (will be only 1 file per pipeline):
// [autosmomal_region, <path_here>/autosomes_N_samples]

ch_N_samples_chrX = Channel.fromFilePairs("${siteqc_results_dir}/*${N_samples_suffix}{_XX,_XY}")
// Must always look as following (will only be one tuple of two files per pipeline):
// [sex_chr, [<path_here>/sex_chr_N_samples_XX, <path_here>/sex_chr_N_samples_XY]]

ch_N_samples = ch_N_samples_autosomes.mix(ch_N_samples_chrX)
// .mix() is joining two channels in one that carries number of samples information for both autosomes
// and X chromosome. It will be a channel with only two elements. Each element is a tuple, where
// first item in tuple defines chromosome type - "autosomal_region" or "sex_chr".
// This will be used before end_aggregate_annotation step to get correct N_samples file to each
// bcf tuple based if bcf is comming from autosomal or chrX region.

Channel.fromPath(params.xx_sample_ids, checkIfExists: true)
                    .set { ch_xx_sample_id }


/*
 * STEP 1 - prep_hwe
 */

process prep_hwe {
    publishDir "${params.outdir}/prep_hwe_pop_files", mode: params.publish_dir_mode

    input:
    file(predicted_ancestries) from ch_infer_ancestry
    file(unrelated_list) from ch_unrelated_list
    file(xx_sample_id) from ch_xx_sample_id

    output:
    file("*_unrelated_pop.keep") into ch_unrelated_by_pop_keep_files
    file("*_unrelated_xx.keep") into ch_unrelated_by_pop_keep_files_xx

    script:
    """
    prep_hwe.R --predicted_ancestries='${predicted_ancestries}' \
               --unrelated_list='${unrelated_list}'

    for pop in AFR EUR SAS EAS; do
        grep -Fwf  ${xx_sample_id} \${pop}_unrelated_pop.keep > \${pop}_unrelated_xx.keep
    done
    """
}

// Joining two pop channels by pop value to have pairs of non-sex/xx keep files for given population.
ch_unrelated_by_pop_keep_files.flatten()
                    .map { pop_file -> [pop_file.getName().split('_').first(), pop_file] }
                    .join(
                        ch_unrelated_by_pop_keep_files_xx.flatten()
                            .map { pop_file -> [pop_file.getName().split('_').first(), pop_file] }
                    )
    // Now combining 4 tuples of pop_keep files with all input bcf files to get all possible combinations of both.
    // This will allow to parallelize p_hwe step per bcf per pop_keep file which are all independent.
                    .combine(ch_bcfs_p_hwe)
                    .set { ch_bcfs_pop_tuple }

/*
 * STEP 2 - p_hwe
 */

process p_hwe {
    publishDir "${params.outdir}/p_hwe_plink_files", mode: params.publish_dir_mode

    input:
    tuple val(population), file(unrealted_pop_keep_file), file(unrealted_pop_keep_file_xx),
            val(region), file(bcf), file(index) from ch_bcfs_pop_tuple

    output:
    tuple val(region), file("*.hwe") into ch_p_hwe
    tuple file("*.log"), file("*.nosex") into ch_p_hwe_logs

    script:
    // If region is from chrX, select the female-only xx keep file.
    keep_file = region =~ /chrX/ ? unrealted_pop_keep_file_xx : unrealted_pop_keep_file
    """
    plink --bcf ${bcf} \
          --hardy midp \
          --keep ${keep_file} \
          --double-id \
          --allow-extra-chr \
          --out ${region}_${population}
    """
}


// Joining all files related to one initial bcf file for end_aggregate_annotation step.

// First, .hwe files emmited separately in ch_p_hwe are groupped together in tuples of
// four files that belong to the same input bcf file using .groupTuple() operator.
// Result - a tuple with two items: 1 - region value, 2 - tuple of four hwe files.

// Second, additional item is added to the tuple on the first position, it is defining
// the chromosome type (autosome or chrX) based on if region name has "chrX" substring.
// This chromosome type item is used to join each autosomal/chrX tuple with appropriate
// N_samples file coming from ch_N_samples. .combine() operaator allows to join single
// autosomal N_count file with each autosomal bcf tuple, same for chrX. After combining
// the chromosome type element is removed because not needed anymore.

// Third, the tuple of four hwe files + N_samples is joined with all metrics files from
// SiteQC pipeline that were generated from same input bcf file.
// All groupping is determined by the region value that is firs element of every channel
// and carries the bcf identity information.
ch_joined_files_to_aggregate = ch_p_hwe.groupTuple()
                                .map { region, hwe_files -> [ region =~ /chrX/ ? "sex_chr" : "autosmomal_region" , region, hwe_files] }
                                .combine(ch_N_samples, by: 0)
                                .map { chr_type, region, hwe_files, N_samples_files -> [region, hwe_files, N_samples_files] }
                                .join(ch_startfiles)
                                .join(ch_miss1)
                                .join(ch_miss2)
                                .join(ch_medianCoverageAll)
                                .join(ch_medianNonMiss)
                                .join(ch_medianGQ)
                                .join(ch_hetAll)
                                .join(ch_hetPass)
                                .join(ch_MendErr)
                                .join(ch_AC_counts)

/*
 * Step 3 - end_aggregate_annotation
 */

process end_aggregate_annotation {
    publishDir "${params.outdir}/Annotation_final/", mode: params.publish_dir_mode

    input:
    tuple val(region),
          file(hwe_files),
          file(N_samples),
          file(startfile),
          file(miss1),
          file(miss2),
          file(medianCoverageAll),
          file(medianNonMiss),
          file(medianGQ),
          file(hetAll),
          file(hetPass),
          file(MendErr),
          file(AC_counts) from ch_joined_files_to_aggregate

    output:
    tuple val(region), file("BCFtools_site_metrics_*.txt.gz"), file("BCFtools_site_metrics_*.txt.gz.tbi") into ch_end_aggr_annotation
    file "Summary_stats/*.txt" optional true into ch_summary_stats

    script:
    // The hwe files are not explicitly provided to R script as arguments.
    // Instead R script just expects four hwe files to be present in wordkir.
    // For the chrX case, we don't even provide the pairs of _XX and _XY files,
    // but the basenames instead. Nextflow will stage the files to wdir, R script
    // will handle the names and find the files.
    startfile_basename = startfile[0].toString().replaceAll("_X[XY]","")
    miss1_basename = miss1[0].toString().replaceAll("_X[XY]","")
    miss2_basename = miss2[0].toString().replaceAll("_X[XY]","")
    medianCoverageAll_basename = medianCoverageAll[0].toString().replaceAll("_X[XY]","")
    medianNonMiss_basename = medianNonMiss[0].toString().replaceAll("_X[XY]","")
    medianGQ_basename = medianGQ[0].toString().replaceAll("_X[XY]","")
    hetAll_basenamne = hetAll[0].toString().replaceAll("_X[XY]","")
    hetPass_basenamne = hetPass[0].toString().replaceAll("_X[XY]","")
    MendErr_basename = MendErr[0].toString().replaceAll("_X[XY]","")
    N_samples_basename = N_samples[0].toString().replaceAll("_X[XY]","")
    """
    mkdir -p "Summary_stats"
    if [[ $region == *"chrX"* ]]; then
        sexChromAnnotation.R \
        ${region} \
        ${startfile_basename} \
        ${miss1_basename} \
        ${miss2_basename} \
        ${medianCoverageAll_basename} \
        ${medianNonMiss_basename} \
        ${medianGQ_basename} \
        ${hetAll_basenamne} \
        ${hetPass_basenamne} \
        ${MendErr_basename} \
        '.' \
        ${N_samples_basename}
    else
        annotatePerChunk.R \
        ${region} \
        ${startfile} \
        ${miss1} \
        ${miss2} \
        ${medianCoverageAll} \
        ${medianNonMiss} \
        ${medianGQ} \
        ${hetAll} \
        ${hetPass} \
        ${MendErr} \
        '.' \
        ${N_samples} \
        ${AC_counts} \
        '.' \
        ${params.king} \
        ${params.aggregate_final}
    fi

    bgzip -f BCFtools_site_metrics_${region}.txt
    tabix -f -s1 -b2 -e2 BCFtools_site_metrics_${region}.txt.gz
    """
}



/*
 * Step 4 - Create new bcf header
 */


process make_header {
    publishDir "${params.outdir}/Additional_header/", mode: params.publish_dir_mode

    input:
    file(all_flags) from ch_summary_stats.collect()
    // From main bcf input cahnnel we need only 1 bcf file to get the example header
    file(bcf) from ch_bcfs_make_header.randomSample(1).map{region, bcf, index -> [bcf]}

    output:
    file("additional_header.txt") into ch_additional_header

    script:
    """
    # Construct a file with unique set of all stats ever seen in all bcf stats files.
    tail -n +2 -q ${all_flags} | cut -f1 | sort | uniq  > "all_seen_flags.txt"

    # Get the original header from example bcf file
    bcftools view -h ${bcf} > "orig.hdr"

    # Clean out INFO lines, these wil be repopulated
    grep -v \\#\\#INFO "orig.hdr" > "orig.hdr2" && \
    mv "orig.hdr2" "orig.hdr"

    # Run the R script to construct the header
    headerbuild.R

    """
}


/*
 * Step 5 - final annotation step
 */

process annotate_bcfs {
    publishDir "${params.outdir}/Final_annotated_BCFs/", mode: params.publish_dir_mode

    input:
    tuple val(region), file(bcf), file(index), file(bcf_site_metrics), file(index2) from ch_bcfs_final_annotation.join(ch_end_aggr_annotation)
    file(additional_header) from ch_additional_header

    output:
    tuple file('*.vcf.gz'), file('*.vcf.gz.csi') into ch_final_annotated_vcfs
    script:
    outfile=file(bcf).getSimpleName()+'.vcf.gz'
    """
    bcftools annotate ${bcf} \
    -x FILTER,^INFO/OLD_MULTIALLELIC,^INFO/OLD_CLUMPED \
    -a ${bcf_site_metrics} \
    -h ${additional_header} \
    -c CHROM,POS,REF,ALT,missingness,medianDepthAll,medianDepthNonMiss,medianGQ,completeGTRatio,MendelSite,ABratio,phwe_afr,phwe_eur,phwe_eas,phwe_sas,FILTER,- | \
    bcftools +fill-tags \
    -o ${outfile} \
    -Oz \
    --threads 16 \
    -- -d -t AC,AC_Hom,AC_Het,AC_Hemi,AN

    bcftools index ${outfile}

    """
}





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
