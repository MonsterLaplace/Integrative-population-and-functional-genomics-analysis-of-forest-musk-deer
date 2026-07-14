#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import gzip
import pandas as pd
from collections import defaultdict

def open_maybe_gzip(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")

def load_site_level_table(site_file, keep_classes=None):
    df = pd.read_csv(site_file, sep="\t", dtype=str)
    required = ['CHROM', 'POS', 'REF', 'ALT', 'class']
    for col in required:
        if col not in df.columns:
            raise ValueError(f"Missing required column in site-level table: {col}")

    if keep_classes is not None:
        df = df[df["class"].isin(keep_classes)].copy()

    site_class = {}
    for _, row in df.iterrows():
        key = (row['CHROM'], row['POS'], row['REF'], row['ALT'])
        site_class[key] = row['class']
    return site_class

def parse_gt(sample_field):
    gt = sample_field.split(":")[0]

    if gt in ("./.", ".|."):
        return "missing"
    elif gt in ("0/0", "0|0"):
        return "hom_ref"
    elif gt in ("0/1", "1/0", "0|1", "1|0"):
        return "het"
    elif gt in ("1/1", "1|1"):
        return "hom_alt"
    else:
        # crude handling of multiallelic or complex genotype
        alleles = gt.replace("|", "/").split("/")
        if "." in alleles:
            return "missing"
        alt_count = sum(1 for a in alleles if a != "0")
        if alt_count == 0:
            return "hom_ref"
        elif alt_count == 1:
            return "het"
        elif alt_count >= 2:
            return "hom_alt"
        return "other"

def main():
    parser = argparse.ArgumentParser(description="Calculate per-sample genetic load from site-level functional variants and VCF.")
    parser.add_argument("-s", "--site-table", required=True, help="Site-level functional variant table")
    parser.add_argument("-v", "--vcf", required=True, help="VCF file (gzipped or plain text)")
    parser.add_argument("-o", "--output", required=True, help="Output TSV")
    parser.add_argument("--classes", nargs="+", default=["LoF", "missense_severe", "missense_all"],
                        help="Variant classes to calculate (default: LoF missense_severe missense_all)")
    args = parser.parse_args()

    site_class = load_site_level_table(args.site_table, keep_classes=args.classes)

    # stats[sample][class]
    stats = defaultdict(lambda: defaultdict(lambda: {
        "n_sites_observed": 0,
        "n_missing": 0,
        "n_hom_ref": 0,
        "n_het": 0,
        "n_hom_alt": 0
    }))

    sample_names = []

    with open_maybe_gzip(args.vcf) as f:
        for line in f:
            if line.startswith("##"):
                continue

            if line.startswith("#CHROM"):
                header = line.rstrip("\n").split("\t")
                sample_names = header[9:]
                continue

            parts = line.rstrip("\n").split("\t")
            chrom, pos, ref, alt = parts[0], parts[1], parts[3], parts[4]

            # skip multi-allelic sites for robustness
            if "," in alt:
                continue

            key = (chrom, pos, ref, alt)
            if key not in site_class:
                continue

            cls = site_class[key]
            sample_fields = parts[9:]

            for sname, sfield in zip(sample_names, sample_fields):
                gt_type = parse_gt(sfield)
                stats[sname][cls]["n_sites_observed"] += 1

                if gt_type == "missing":
                    stats[sname][cls]["n_missing"] += 1
                elif gt_type == "hom_ref":
                    stats[sname][cls]["n_hom_ref"] += 1
                elif gt_type == "het":
                    stats[sname][cls]["n_het"] += 1
                elif gt_type == "hom_alt":
                    stats[sname][cls]["n_hom_alt"] += 1

    rows = []
    for sample in sorted(stats.keys()):
        for cls in sorted(stats[sample].keys()):
            rec = stats[sample][cls]
            denominator = 2 * rec["n_hom_alt"] + rec["n_het"]
            if denominator > 0:
                genetic_load = 2 * rec["n_hom_alt"] / denominator
            else:
                genetic_load = None

            rows.append({
                "sample": sample,
                "class": cls,
                "n_sites_observed": rec["n_sites_observed"],
                "n_missing": rec["n_missing"],
                "n_hom_ref": rec["n_hom_ref"],
                "n_het": rec["n_het"],
                "n_hom_alt": rec["n_hom_alt"],
                "genetic_load": genetic_load
            })

    out = pd.DataFrame(rows)
    out.to_csv(args.output, sep="\t", index=False)

    print("=== Genetic load summary written ===")
    print(args.output)

if __name__ == "__main__":
    main()
