#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pandas as pd
import numpy as np

def load_fai(fai_file):
    fai = pd.read_csv(
        fai_file,
        sep="\t",
        header=None,
        usecols=[0, 1],
        names=["CHROM", "length"]
    )
    return fai

def make_windows(fai, window_size=50000, step_size=10000):
    rows = []
    for _, r in fai.iterrows():
        chrom = str(r["CHROM"])
        chrom_len = int(r["length"])

        for start in range(1, chrom_len + 1, step_size):
            end = start + window_size - 1
            if end > chrom_len:
                end = chrom_len
            rows.append({
                "CHROM": chrom,
                "start": start,
                "end": end
            })
    return pd.DataFrame(rows)

def load_plink_hom(hom_file):
    # 用 sep=r"\s+" 代替 delim_whitespace=True，兼容更多 pandas 版本
    df = pd.read_csv(hom_file, sep=r"\s+", engine="python", dtype=str)

    required_cols = ["IID", "CHR", "POS1", "POS2"]
    for c in required_cols:
        if c not in df.columns:
            raise ValueError(f"Missing required column in .hom file: {c}")

    df["CHR"] = df["CHR"].astype(str)
    df["POS1"] = pd.to_numeric(df["POS1"], errors="coerce")
    df["POS2"] = pd.to_numeric(df["POS2"], errors="coerce")
    df = df.dropna(subset=["POS1", "POS2"])

    df["POS1"] = df["POS1"].astype(int)
    df["POS2"] = df["POS2"].astype(int)

    # 保证 roh_start <= roh_end
    df["roh_start"] = df[["POS1", "POS2"]].min(axis=1)
    df["roh_end"] = df[["POS1", "POS2"]].max(axis=1)

    return df[["IID", "CHR", "roh_start", "roh_end"]].copy()

def compute_roh_window_frequency(windows, roh_df, total_samples):
    results = []

    # 按染色体拆分
    roh_by_chr = {chrom: sub.copy() for chrom, sub in roh_df.groupby("CHR")}
    win_by_chr = {chrom: sub.copy() for chrom, sub in windows.groupby("CHROM")}

    for chrom, wsub in win_by_chr.items():
        rsub = roh_by_chr.get(chrom, pd.DataFrame(columns=["IID", "CHR", "roh_start", "roh_end"]))

        if rsub.empty:
            for _, w in wsub.iterrows():
                results.append({
                    "CHROM": w["CHROM"],
                    "start": int(w["start"]),
                    "end": int(w["end"]),
                    "roh_freq": 0.0,
                    "mean_roh_count": 0.0
                })
            continue

        for _, w in wsub.iterrows():
            w_start = int(w["start"])
            w_end = int(w["end"])

            # 找与当前窗口重叠的 ROH
            overlaps = rsub[(rsub["roh_end"] >= w_start) & (rsub["roh_start"] <= w_end)]

            n_samples_hit = overlaps["IID"].nunique()
            total_roh_count = overlaps.shape[0]

            roh_freq = n_samples_hit / total_samples if total_samples > 0 else np.nan
            mean_roh_count = total_roh_count / total_samples if total_samples > 0 else np.nan

            results.append({
                "CHROM": chrom,
                "start": w_start,
                "end": w_end,
                "roh_freq": roh_freq,
                "mean_roh_count": mean_roh_count
            })

    out = pd.DataFrame(results)
    return out

def main():
    parser = argparse.ArgumentParser(
        description="Generate roh_window_frequency.tsv from PLINK .hom and reference .fai"
    )
    parser.add_argument("-r", "--roh-hom", required=True, help="PLINK .hom file")
    parser.add_argument("-f", "--fai", required=True, help="reference genome .fai file")
    parser.add_argument("-o", "--output", required=True, help="output roh_window_frequency.tsv")
    parser.add_argument("--window", type=int, default=50000, help="window size (default 50000)")
    parser.add_argument("--step", type=int, default=10000, help="step size (default 10000)")
    args = parser.parse_args()

    fai = load_fai(args.fai)
    windows = make_windows(fai, window_size=args.window, step_size=args.step)
    roh_df = load_plink_hom(args.roh_hom)

    total_samples = roh_df["IID"].nunique()
    out = compute_roh_window_frequency(windows, roh_df, total_samples)

    out.to_csv(args.output, sep="\t", index=False)

    print("=== ROH window frequency summary ===")
    print(f"Total windows: {out.shape[0]}")
    print(f"Total samples with at least one ROH: {total_samples}")
    print(f"Windows with roh_freq > 0: {(out['roh_freq'] > 0).sum()}")
    print(f"Windows with roh_freq >= 0.5: {(out['roh_freq'] >= 0.5).sum()}")
    print(f"Output written to: {args.output}")

if __name__ == "__main__":
    main()
