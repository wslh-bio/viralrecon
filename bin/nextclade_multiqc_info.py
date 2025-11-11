#!/usr/bin/env python3
"""
Edit Nextclade TSV file for MultiQC to include the database name,
version, and reference accession number for each sample processed by nextclade_run.

Usage:
    nextclade_multiqc_info.py --db_dir <db_dir> --clade_tsv <clade_tsv> --out_tsv <out_tsv>

Reads the Nextclade TSV file and adds database version information from pathogen.json.
"""

import argparse
import json
import logging
import os
import pandas as pd
import shutil
import sys
import zipfile

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")


def prepare_db_dir(db_dir):
    """Extract ZIP if provided, else return directory path."""
    if zipfile.is_zipfile(db_dir):
        folder_name = os.path.splitext(os.path.basename(db_dir))[0]
        shutil.rmtree(folder_name, ignore_errors=True)
        logging.info("Extracting database ZIP %s to %s", db_dir, folder_name)
        with zipfile.ZipFile(db_dir, "r") as zip_ref:
            zip_ref.extractall(folder_name)
        return folder_name
    return db_dir


def find_pathogen_json(db_path):
    """Locate pathogen.json in the given directory."""
    json_file = os.path.join(db_path, "pathogen.json")
    if os.path.isfile(json_file):
        return json_file
    for root, _, files in os.walk(db_path):
        if "pathogen.json" in files:
            return os.path.join(root, "pathogen.json")
    return None


def build_db_info(json_path):
    """Extract dataset info from pathogen.json."""
    with open(json_path, encoding="utf-8") as fh:
        data = json.load(fh)

    name = data.get("attributes", {}).get("name") or data.get("name", "Unknown")
    ref_name = data.get("attributes", {}).get("reference name") or data.get("reference", {}).get("name", "Unknown")
    ref_acc = data.get("attributes", {}).get("reference accession") or data.get("reference", {}).get(
        "accession", "Unknown"
    )

    version = None
    if isinstance(data.get("version"), dict):
        version = data.get("version").get("tag") or data.get("version").get("name")
    else:
        version = data.get("version")

    return {
        "nextclade_dataset": {
            "name": name,
            "reference": {"name": ref_name, "accession": ref_acc},
            "version": version or "Unknown",
        }
    }


def update_tsv(clade_tsv, out_tsv, db_info):
    """Append database info columns to TSV."""
    if not os.path.isfile(clade_tsv):
        logging.error("Input TSV file %s not found", clade_tsv)
        return 1

    logging.info("Reading TSV with pandas")
    try:
        df = pd.read_csv(clade_tsv, sep="\t")
    except Exception as e:
        logging.error("Failed to read TSV: %s", e)
        return 6

    if df.empty or len(df.columns) < 2:
        logging.error("Invalid TSV: empty or missing header")
        return 5

    db = db_info.get("nextclade_dataset", {})
    logging.info("Adding database info columns")
    df["nextclade_db_name"] = db.get("name", "Unknown")
    df["nextclade_db_version"] = db.get("version", "Unknown")
    df["nextclade_db_reference_name"] = db.get("reference", {}).get("name", "Unknown")
    df["nextclade_db_reference_accession"] = db.get("reference", {}).get("accession", "Unknown")

    logging.info("Writing updated TSV to %s", out_tsv)
    try:
        df.to_csv(out_tsv, sep="\t", index=False)
    except Exception as e:
        logging.error("Failed to write TSV: %s", e)
        return 7

    logging.info("TSV update completed successfully.")
    return 0


class CompileResults(argparse.ArgumentParser):
    """Custom ArgumentParser that prints help on error."""

    def error(self, message):
        self.print_help()
        sys.stderr.write(f"\nERROR DETECTED: {message}\n")
        sys.exit(1)


if __name__ == "__main__":
    parser = CompileResults(
        prog="nextclade_multiqc_info.py",
        description="Update Nextclade TSV with database info for MultiQC",
        epilog="Example usage: python nextclade_multiqc_info.py --db_dir <db_dir> --clade_tsv <clade_tsv> --out_tsv <out_tsv>",
    )

    parser.add_argument("--db_dir", required=True, help="Directory containing pathogen.json or ZIP")
    parser.add_argument("--clade_tsv", required=True, help="Input Nextclade TSV file")
    parser.add_argument("--out_tsv", required=True, help="Output TSV file")

    args = parser.parse_args()

    logging.info("Preparing database directory")
    db_dir = prepare_db_dir(args.db_dir)

    logging.info("Finding pathogen.json")
    json_file = find_pathogen_json(db_dir)
    if json_file is None:
        logging.error("pathogen.json not found in %s", db_dir)
        sys.exit(2)

    logging.info("Building database info")
    info = build_db_info(json_file)

    logging.info("Updating TSV file with pandas")
    sys.exit(update_tsv(args.clade_tsv, args.out_tsv, info))
