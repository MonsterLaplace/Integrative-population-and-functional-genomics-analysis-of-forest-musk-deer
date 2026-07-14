#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import gzip
import re
import pandas as pd
from collections import defaultdict

def open_maybe_gzip(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")

# --------------------------------------------------
# Chromosome name normalization
# Goal: convert all styles to unified form:
#   Chr01, Chr02, ..., Chr10, ..., ChrX
# Accept examples:
#   1 -> Chr01
#   2 -> Chr02
#   X -> ChrX
#   chr1 -> Chr01
#   Chr1 -> Chr01
#   Chr01 -> Chr01
#   chr01 -> Chr01
#   ChrX -> ChrX
# --------------------------------------------------
def normalize_chrom(chrom):
    if chrom is None:
        return None

    c = str(chrom).strip()
    if c == "":
        return None

    # remove common prefixes
    c = re.sub(r'^(chr|Chr)', '', c)

    # if purely numeric
    if re.fullmatch(r'\d+', c):
        num = int(c)
        return f"Chr{num:02d}"

    # if sex chromosome X/Y or mt-like
    cu = c.upper()
    if cu == "X":
        return "ChrX"
    if cu == "Y":
        return "ChrY"
    if cu in ["M", "MT", "MITO", "MITOCHONDRIA"]:
        return "ChrM"

    # fallback: preserve original label style after removing prefix
    return "Chr" + c

def load_site_classes(site_file, keep_classes=None):
    df = pd.read_csv(site_file, sep="\t", dtype=str)

    required = ["CHROM", "POS", "REF", "ALT", "class"]
    for col in required:
        if col not in df.columns:
            raise ValueError(f"Missing required column in site-level table: {col}")

    if keep_classes is not None:
        df = df[df["class"].isin(keep_classes)].copy()

    site_class = {}
    for _, row in df.iterrows():
        chrom = normalize_chrom(row["CHROM"])
        key = (chrom, int(row["POS"]), row["REF"], row["ALT"])
        site_class[key] = row["class"]

    return site_class

def load_roh_bed(roh_bed):
    """
    ROH BED format:
    chr   start   end   sample
    BED start is 0-based, end is 1-based open interval
    convert to 1-based inclusive coordinates for VCF comparison
    """
    roh_by_sample_chr = defaultdict(lambda: defaultdict(list))

    with open(roh_bed) as f:
        for line_num, line in enumerate(f, start=1):
            if not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) < 4:
                raise ValueError(f"ROH BED line {line_num} has fewer than 4 columns")

            chrom_raw, start, end, sample = parts[:4]
            chrom = normalize_chrom(chrom_raw)
            start = int(start) + 1   # BED -> 1-based inclusive
            end = int(end)

            roh_by_sample_chr[sample][chrom].append((start, end))

    return roh_by_sample_chr

def in_any_roh(pos, roh_list):
    for s, e in roh_list:
        if s <= pos <= e:
            return True
    return False

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
        # crude handling for multiallelic/complex genotypes
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
    parser = argparse.ArgumentParser(
        description="Calculate genetic load inside vs outside ROH, with chromosome name normalization."
    )
    parser.add_argument("-s", "--site-table", required=True, help="Site-level functional variant table")
    parser.add_argument("-v", "--vcf", required=True, help="VCF file (plain or gzipped)")
    parser.add_argument("-r", "--roh-bed", required=True, help="ROH BED file: chr start end sample")
    parser.add_argument("-o", "--output", required=True, help="Output TSV")
    parser.add_argument("--classes", nargs="+",
                        default=["LoF", "missense_severe", "missense_all"],
                        help="Variant classes to include (default: LoF missense_severe missense_all)")
    args = parser.parse_args()

    print("[INFO] Loading site-level variant classes...")
    site_class = load_site_classes(args.site_table, keep_classes=args.classes)

    print("[INFO] Loading ROH BED...")
    roh_by_sample_chr = load_roh_bed(args.roh_bed)

    # stats[sample][class][region]
    stats = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: {
        "n_sites_observed": 0,
        "n_missing": 0,
        "n_hom_ref": 0,
        "n_het": 0,
        "n_hom_alt": 0
    })))

    sample_names = []
    matched_sites = 0
    total_variant_lines = 0

    print("[INFO] Parsing VCF...")
    with open_maybe_gzip(args.vcf) as f:
        for line in f:
            if line.startswith("##"):
                continue

            if line.startswith("#CHROM"):
                header = line.rstrip("\n").split("\t")
                sample_names = header[9:]
                print(f"[INFO] Found {len(sample_names)} samples in VCF")
                continue

            parts = line.rstrip("\n").split("\t")
            chrom_raw, pos, ref, alt = parts[0], int(parts[1]), parts[3], parts[4]
            total_variant_lines += 1

            # skip multi-allelic sites
            if "," in alt:
                continue

            chrom = normalize_chrom(chrom_raw)
            key = (chrom, pos, ref, alt)
            if key not in site_class:
                continue

            matched_sites += 1
            cls = site_class[key]
            sample_fields = parts[9:]

            for sname, sfield in zip(sample_names, sample_fields):
                roh_list = roh_by_sample_chr.get(sname, {}).get(chrom, [])
                region = "inside_ROH" if in_any_roh(pos, roh_list) else "outside_ROH"

                gt_type = parse_gt(sfield)
                stats[sname][cls][region]["n_sites_observed"] += 1

                if gt_type == "missing":
                    stats[sname][cls][region]["n_missing"] += 1
                elif gt_type == "hom_ref":
                    stats[sname][cls][region]["n_hom_ref"] += 1
                elif gt_type == "het":
                    stats[sname][cls][region]["n_het"] += 1
                elif gt_type == "hom_alt":
                    stats[sname][cls][region]["n_hom_alt"] += 1

    rows = []
    for sample in sorted(stats.keys()):
        for cls in sorted(stats[sample].keys()):
            # 保证 inside/outside 都输出
            for region in ["inside_ROH", "outside_ROH"]:
                rec = stats[sample][cls][region]
                denominator = 2 * rec["n_hom_alt"] + rec["n_het"]

                if denominator > 0:
                    genetic_load = 2 * rec["n_hom_alt"] / denominator
                else:
                    genetic_load = None

                rows.append({
                    "sample": sample,
                    "class": cls,
                    "region": region,
                    "n_sites_observed": rec["n_sites_observed"],
                    "n_missing": rec["n_missing"],
                    "n_hom_ref": rec["n_hom_ref"],
                    "n_het": rec["n_het"],
                    "n_hom_alt": rec["n_hom_alt"],
                    "genetic_load": genetic_load
                })

    out = pd.DataFrame(rows)
    out.to_csv(args.output, sep="\t", index=False)

    print(f"[INFO] Total VCF variant lines scanned: {total_variant_lines}")
    print(f"[INFO] Matched functional sites: {matched_sites}")
    print(f"[INFO] Output written to: {args.output}")

if __name__ == "__main__":
    main()
