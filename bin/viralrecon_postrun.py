#!/usr/bin/env python

# Import Dependencies
import glob
import os
import sys
import argparse
import pandas as pd
import pysam
from os.path import exists

# Parse Arguements
def parse_args(args=None):
    """Parse args function"""
    description = (
        "This script creates a summary output file containing other informative output metrics appended to the summary variants metrics file generated by nf-core viralrecon."
    )
    epilog = (
        "Example usasge: python viralrecon_postrun.py [args...]"
    )
    parser = argparse.ArgumentParser(description=description, epilog=epilog)
    parser.add_argument(
        "-s",
        "--summary",
        type=str,
        dest="SUMMARY",
        default="",
        help="The path to variants metric summary.",
    )
    parser.add_argument(
        "-b",
        "--bam_files",
        type=str,
        nargs='+',
        dest="BAM_FILES",
        default="",
        help="Paths to bam files.",
    )
    parser.add_argument(
        "-p",
        "--panoglin",
        type=str,
        nargs='+',
        dest="PANGOLIN",
        default="",
        help="The path to pangolin reports.",
    )
    parser.add_argument(
        "-n",
        "--nextclade",
        type=str,
        nargs='+',
        dest="NEXTCLADE",
        default="",
        help="The path to the nextclade reports.",
    )
    parser.add_argument(
        "-wv",
        "--workflow_version",
        type=str,
        nargs='+',
        action='store',
        dest="WORKFLOW_VERSION",
        default="",
        help="The version of the workflow.",
    )
    parser.add_argument(
        "-rn",
        "--run_name",
        type=str,
        nargs='+',
        action='store',
        dest="RUN_NAME",
        default="",
        help="The name of the run.",
    )
    return parser.parse_args(args)

def get_summary_file(summary):
    sum_df = pd.read_csv(summary)
    renamed = "wslh_" + summary
    sum_df.to_csv(renamed, sep=',', index=False)
    return renamed


def get_pangolin_data(pangolin_list):
    # takes all of the pangolin.csv files into a list and puts them into a combined dataframe
    df_list = []
    for file in pangolin_list:
        item_df = pd.read_csv(file)
        df_list.append(item_df)
    pangolin_df = pd.concat(df_list)
    pangolin_df.to_csv("int_combined_pangolin.csv", index=False)

    """ in order to generate the final postrun output, we need to be able to join dataframes on sample name. To do this, we need to
    parse out the reference fasta file name from the sample name of the pangolin output.
    """

    with open("int_combined_pangolin.csv", 'r') as io, open("final_combined_pangolin.csv", 'w') as fo:
        line_dict = {}
        for line in io:
            if line.startswith("taxon"):  # keep the header line for our final combined pangolin output file
                header = line
            else:
                p_sample = line.split(",")[0]
                sample = p_sample.split(" ")[0].split("/")[0]  # second split is to cover nanopore workflow, since it adds forwards slashes to sample names
                line_dict[sample] = ",".join(line.split(",")[1:])  # the dict has desired sample ID as the key and rest of line as value
        fo.write(header)
        for k, v in line_dict.items():
            fo.write(k + "," + v)
    return "final_combined_pangolin.csv"


def samtools_coverage(bam_list):
    cov_file_list = []
    cov_df_list = []
    # run samtools coverage on each BAM file. This gets "depth after trimming" and "1X coverage after trimming"
    for bam in bam_list:
        samplename = bam.split(".bam")[0].rstrip()
        coverage_file = samplename + ".cov.txt"
        try:
            pysam.coverage(bam, "-o", coverage_file)
            cov_file_list.append(coverage_file)
        except:
            pass

    try:
        for c_file in cov_file_list:
                sample_id = c_file.split(".")[0]
                cfile_df = pd.read_csv(c_file, sep="\t")
                cfile_df.insert(0, "Sample", sample_id)
                cov_df_list.append(cfile_df)
                coverage_df = pd.concat(cov_df_list)
        coverage_df.to_csv("combined_samtools_cov_report.csv", index=False)
    except:
        pass
    return "combined_samtools_cov_report.csv"


def get_nextclade_data(nextclade_list):
    # same strategy as 'pangolin' above
    df_list = []
    for file in nextclade_list:
        item_df = pd.read_csv(file, sep =";")
        df_list.append(item_df)
    nextclade_df = pd.concat(df_list)
    if 'index' in nextclade_df.columns:
        nextclade_df = nextclade_df.drop(columns=['index'])
    nextclade_df.to_csv("int_combined_nextclade.csv", index=False)


    with open("int_combined_nextclade.csv", 'r') as io, open("final_combined_nextclade.csv", 'w') as fo:
        line_dict = {}
        for line in io:
            if line.startswith("seqName"):
                header = line
            else:
                sample = line.split(" ")[0].split("/")[0]  # second split is to cover nanopore workflow, since it adds forwards slashes to sample names
                line_dict[sample] = ",".join(line.split(",")[1:]) 
        fo.write(header)
        for k, v in line_dict.items():
            fo.write(k + "," + v)
    return "final_combined_nextclade.csv"

def combine_results(summary_file, samtools_file, pangolin_file, nextclade_file, workflow_version, wslh_output):
    # function that takes in any files generated from this script and adds to the summary dataframe
    sum_df = pd.read_csv(summary_file, dtype = str)
    sum_df['depth_after_trimming'] = None
    sum_df['1X_coverage_after_trimming'] = None
    sum_df['viralrecon_version'] = None
    columns = list(sum_df.columns)
    columns.remove('Sample')


    if exists(samtools_file):
        samtools_df = pd.read_csv(samtools_file, dtype=str, usecols=['Sample', 'coverage', 'meandepth'], index_col=False)
        samtools_df = samtools_df.add_prefix('samtools_')
        st_cols = list(samtools_df.columns)
        st_cols.remove('samtools_Sample')
        st_cols.remove('samtools_meandepth')
        st_cols.remove('samtools_coverage')


        sum_df = pd.merge(sum_df, samtools_df, left_on = "Sample", right_on="samtools_Sample", how = 'outer')
        sum_df['depth_after_trimming'].fillna(sum_df['samtools_meandepth'], inplace=True)
        sum_df['1X_coverage_after_trimming'].fillna(sum_df['samtools_coverage'], inplace=True)
        columns = ['Sample'] + columns
    else:
        pass

    if exists(pangolin_file):
        p_df = pd.read_csv(pangolin_file, dtype=str)
        p_df = p_df.add_prefix('pangolin_')
        p_cols = list(p_df.columns)
        p_cols.remove('pangolin_taxon')
        p_cols.remove('pangolin_lineage')

        sum_df = pd.merge(sum_df, p_df, left_on="Sample", right_on="pangolin_taxon", how = 'outer')
        sum_df.drop('pangolin_taxon', axis=1, inplace=True)
        columns = columns + p_cols
    else:
        pass

    if exists(nextclade_file):
        n_df = pd.read_csv(nextclade_file, dtype=str, usecols=['seqName', 'qc.overallScore', 'qc.overallStatus'])
        n_df = n_df.add_prefix('nextclade_')
        n_cols = list(n_df.columns)
        n_cols.remove('nextclade_seqName')

        sum_df = pd.merge(sum_df, n_df, left_on="Sample", right_on="nextclade_seqName", how = 'outer')
        columns = columns + n_cols
    else:
        pass

    sum_df["sample"] = sum_df["Sample"].str.split('_S').str[0]
    sum_df = sum_df.assign(viralrecon_version=workflow_version[0])
    df = sum_df[['sample', 'Sample', 'Pangolin lineage', 'Nextclade clade', 'pangolin_scorpio_call', 'pangolin_qc_status',
                 '# Input reads', '# Trimmed reads (fastp)', '% Non-host reads (Kraken 2)',
                 '% Mapped reads', '# Mapped reads', 'depth_after_trimming', '1X_coverage_after_trimming',
                 'Coverage median', '% Coverage > 1x', '% Coverage > 10x', '# SNPs', '# INDELs', '# Missense variants',
                 '# Ns per 100kb consensus', '# Trimmed reads (iVar)', 'nextclade_qc.overallStatus', 'pangolin_conflict', 'pangolin_ambiguity_score',
                 'pangolin_scorpio_support', 'pangolin_scorpio_conflict', 'pangolin_scorpio_notes', 'pangolin_version',
                 'pangolin_pangolin_version', 'pangolin_scorpio_version', 'pangolin_constellation_version',
                 'pangolin_is_designated', 'nextclade_qc.overallScore', 'pangolin_qc_notes', 'pangolin_note', 'viralrecon_version']]
    df1 = df.rename(columns={'sample':'sample_id', 'Sample':'sample', 'Nextclade clade':'nextclade_clade', 'Pangolin lineage':'pangolin_lineage',
                       '# Input reads':'num_raw_reads', '# Trimmed reads (fastp)': 'num_trimmed_reads_fastp',
                       '# Trimmed reads (iVar)':'num_trimmed_reads_iVar', '% Non-host reads (Kraken 2)':'percent_non-host_reads_kraken2',
                       '% Mapped reads':'percent_mapped_reads', '# Mapped reads':'num_mapped_reads',
                       'Coverage median':'coverage_median', '% Coverage > 1x':'percent_coverage>1x',
                       '% Coverage > 10x':'percent_coverage>10x', '# SNPs':'num_SNPs', '# INDELs':'num_INDELs',
                       '# Missense variants':'num_missense_variants', '# Ns per 100kb consensus':'num_Ns_per_100kb_consensus' })
    df1.to_csv(wslh_output[0]+"_viralrecon_report.csv", index=False)

    return wslh_output[0]+"_viralrecon_report.csv"

def wslh_qc_filter(viralrecon_report, wslh_output):
    # Read the spreadsheet into a DataFrame
    df = pd.read_csv(viralrecon_report)

    # Apply criteria to fill the 'WSLH_qc' column
    mask_ntc = (
        ((df['sample_id'].str.contains('VQ')) & (df['1X_coverage_after_trimming'] < 10)) |
        ((df['1X_coverage_after_trimming'].isna()))
    )

    mask_other = (
        (~mask_ntc) &
        (df['depth_after_trimming'] >= 100) &
        (df['1X_coverage_after_trimming'] >= 90) &
        (df['pangolin_lineage'].notna()) &
        (df['pangolin_scorpio_call'].notna())
    )
    
    df.loc[(mask_ntc & df['sample_id'].str.contains('VQ')), 'WSLH_qc'] = 'NTC_pass'
    df.loc[(~mask_ntc & df['sample_id'].str.contains('VQ')), 'WSLH_qc'] = 'NTC_fail'
    df.loc[(mask_other & df['sample_id'].str.contains('VR')), 'WSLH_qc'] = 'pass'
    df.loc[(~mask_other & df['sample_id'].str.contains('VR')), 'WSLH_qc'] = 'fail'

    # Set 'unassigned_failed_qc' for rows with blanks in 'pangolin_lineage'
    #df.loc[df['pangolin_lineage'].isna() & df['sample_id'].str.contains('VR'), 'pangolin_lineage'] = 'unassigned_failed_qc'

    # Reindex columns
    df = df.reindex(columns = ['sample_id', 'WSLH_qc', 'pangolin_lineage', 'nextclade_clade', 'pangolin_scorpio_call', 'pangolin_qc_status', 'num_raw_reads', 'num_trimmed_reads_fastp', 'percent_non-host_reads_kraken2', 'percent_mapped_reads', 'num_mapped_reads', 'depth_after_trimming', '1X_coverage_after_trimming', 'coverage_median', 'percent_coverage>1x', 'percent_coverage>10x', 'num_SNPs', 'num_INDELs', 'num_missense_variants', 'num_Ns_per_100kb_consensus', 'num_trimmed_reads_iVar', 'nextclade_qc.overallStatus', 'pangolin_conflict', 'pangolin_ambiguity_score', 'pangolin_scorpio_support', 'pangolin_scorpio_conflict', 'pangolin_scorpio_notes', 'pangolin_version', 'pangolin_pangolin_version', 'pangolin_scorpio_version', 'pangolin_constellation_version', 'pangolin_is_designated', 'nextclade_qc.overallScore', 'pangolin_qc_notes', 'pangolin_note', 'viralrecon_version', 'sample'])

    # Save the updated DataFrame to a new spreadsheet
    df.to_csv(wslh_output[0]+"_viralrecon_report.csv", index=False)

def main(args=None):
    args = parse_args(args)

    sum_file = get_summary_file(summary=args.SUMMARY)
    pangolin_file = get_pangolin_data(pangolin_list=args.PANGOLIN)
    nextclade_file = get_nextclade_data(nextclade_list=args.NEXTCLADE)
    samtools_file = samtools_coverage(bam_list=args.BAM_FILES)
    report = combine_results(sum_file, samtools_file, pangolin_file, nextclade_file, workflow_version=args.WORKFLOW_VERSION, wslh_output=args.RUN_NAME)
    wslh_qc_filter(report, wslh_output=args.RUN_NAME)
if __name__ == "__main__":
    sys.exit(main())