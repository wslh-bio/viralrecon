#!/usr/bin/env python3
"""
Edit Nextclade TSV file for MultiQC to include the database name,
version, and reference accession number for each sample processed by nextclade_run.

Usage:
  nextclade_multiqc_info.py --db_dir <db_dir> --clade_tsv <clade_tsv> --out_tsv <out_tsv>

Reads the Nextclade TSV file and adds database version information from pathogen.json.
"""

import os
import sys
import json
import argparse
import zipfile
import shutil


def parse_args():
    parser = argparse.ArgumentParser(
        description="Augment Nextclade TSV with database info for MultiQC"
    )
    parser.add_argument("--db_dir", required=True, help="Directory containing pathogen.json")
    parser.add_argument("--clade_tsv", required=True, help="Input Nextclade TSV file")
    parser.add_argument("--out_tsv", required=True, help="Output TSV file")
    return parser.parse_args()



def prepare_db_dir(db_dir):
    if zipfile.is_zipfile(db_dir):
        folder_name = os.path.splitext(os.path.basename(db_dir))[0]
        shutil.rmtree(folder_name, ignore_errors=True)
        with zipfile.ZipFile(db_dir, 'r') as zip_ref:
            zip_ref.extractall(folder_name)
        return folder_name
    else:
        return db_dir

def find_pathogen_json(db_path):
    json_file = os.path.join(db_path, "pathogen.json")
    if os.path.isfile(json_file):
        return json_file
    for root, _, files in os.walk(db_path):
        if "pathogen.json" in files:
            return os.path.join(root, "pathogen.json")
    return None


def build_db_info(json_path):
    with open(json_path) as fh:
        data = json.load(fh)
    name = data.get("attributes", {}).get("name") or data.get("name")
    ref_name = (
        data.get("attributes", {}).get("reference name")
        or (data.get("reference", {}).get("name") if data.get("reference") else None)
    )
    ref_acc = (
        data.get("attributes", {}).get("reference accession")
        or (data.get("reference", {}).get("accession") if data.get("reference") else None)
    )
    version = None
    if isinstance(data.get("version"), dict):
        version = data.get("version").get("tag") or data.get("version").get("name")
    else:
        version = data.get("version")
    return {
        "nextclade_dataset": {
            "name": name or "Unknown",
            "reference": {
                "name": ref_name or "Unknown",
                "accession": ref_acc or "Unknown",
            },
            "version": version or "Unknown",
        }
    }


def update_tsv(clade_tsv, out_tsv, db_info):
    db = db_info.get("nextclade_dataset", {})
    dataset_name = db.get("name", "Unknown")
    version = db.get("version", "Unknown")
    reference = db.get("reference", {})
    reference_name = reference.get("name", "Unknown")
    reference_accession = reference.get("accession", "Unknown")

    db_columns = [
        "nextclade_db_name",
        "nextclade_db_version",
        "nextclade_db_reference_name",
        "nextclade_db_reference_accession",
    ]
    db_values = [dataset_name, version, reference_name, reference_accession]

    if not os.path.isfile(clade_tsv):
        sys.stderr.write(f"Error: Input TSV file {clade_tsv} not found\n")
        return 1

    with open(clade_tsv) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        rows = [line.rstrip("\n").split("\t") for line in fh if line.strip()]

    with open(out_tsv, "w") as fh:
        fh.write("\t".join(header + db_columns) + "\n")
        for row in rows:
            if len(row) < len(header):
                row += [""] * (len(header) - len(row))
            fh.write("\t".join(row + db_values) + "\n")

    return 0


def main():
    args = parse_args()
    db_dir = prepare_db_dir(args.db_dir)
    clade_tsv = args.clade_tsv
    out_tsv = args.out_tsv

    json_file = find_pathogen_json(db_dir)
    if json_file is None:
        sys.stderr.write(f"Error: pathogen.json not found in {db_dir}\n")
        return 2

    info = build_db_info(json_file)
    return update_tsv(clade_tsv, out_tsv, info)


if __name__ == "__main__":
    sys.exit(main())
