#!/usr/bin/env python3

import argparse
import logging
import sys
import pandas as pd

logging.basicConfig(level=logging.INFO, format='%(levelname)s : %(message)s')

def add_run_metadata(input_csv, run_name, wf_version):
    logging.debug("Adding run metadata to dataframe")
    df = pd.read_csv(input_csv, sep=',')

    df['run_name'] = run_name
    df['workflow_version'] = wf_version
    
    output_csv = f"{run_name}_viralrecon_report.csv"
    df.to_csv(output_csv, index=False)
    logging.info(f"Wrote {output_csv} with run metadata added to {input_csv}")

class CompileResults(argparse.ArgumentParser):
    def error(self, message):
        self.print_help()
        sys.stderr.write(f'\nERROR DETECTED: {message}\n')
        sys.exit(1)

if __name__ == "__main__":
    parser = CompileResults(
        prog="Compiles WSLH report",
        description="Add run and workflow metadata columns to a CSV report.",
        epilog = "Example usage: python wslh_report_summary.py -s <csv_variants> -r <RUNNAME> -wv <WORKFLOW_VERSION>"
        )
    parser.add_argument(
        "-s",
        "--summary_report",
        help="File containing QC stats in csv format."
    )
    parser.add_argument(
        "-r",
        "--run_name",
        type=str,
        help="Run name or batch of Candida auris in format CA_<machine>_YYMMDD.",
    )
    parser.add_argument(
        "-wv",
        "--wf-version",
        type=str,
        help="Version of the workflow used to generate the results.",
    )

    logging.debug("Run parser to call arguments downstream")
    args = parser.parse_args()

    logging.info("Create wslh report from input files")
    add_run_metadata(args.summary_report, args.run_name, args.wf_version)
