#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pandas as pd


def find_column(df, candidates):
    """在DataFrame列名中查找候选列，忽略大小写"""
    lower_map = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in lower_map:
            return lower_map[cand.lower()]
    return None


def main():
    parser = argparse.ArgumentParser(description="Aggregate iHS SNP scores to windows.")
    parser.add_argument("-i", "--input", required=True, help="normalized iHS file")
    parser.add_argument("-o", "--output", required=True, help="output window table")
    parser.add_argument("--window", type=int, default=50000, help="window size")
    parser.add_argument("--step", type=int, default=10000, help="step size")
    args = parser.parse_args()

    # 读取带表头文件
    df = pd.read_csv(args.input, sep=r"\s+", engine="python")

    # 自动识别列
    chrom_col = find_column(df, ["chr", "chrom", "chromosome"])
    snp_col   = find_column(df, ["id", "snp", "marker"])
    pos_col   = find_column(df, ["pos", "position", "bp"])
    score_col = find_column(df, ["norm_ihs", "ihs"])
    crit_col  = find_column(df, ["crit"])

    if chrom_col is None or pos_col is None or score_col is None:
        raise ValueError(
            f"Unexpected iHS file format. Need columns like chr, pos, norm_ihs/ihs.\n"
            f"Detected columns: {list(df.columns)}"
        )

    # 标准化列名
    keep_cols = [chrom_col, pos_col, score_col]
    rename_map = {
        chrom_col: "CHROM",
        pos_col: "POS",
        score_col: "SCORE"
    }

    if snp_col is not None:
        keep_cols.append(snp_col)
        rename_map[snp_col] = "SNP"
    if crit_col is not None:
        keep_cols.append(crit_col)
        rename_map[crit_col] = "CRIT"

    df = df[keep_cols].rename(columns=rename_map)

    df["POS"] = pd.to_numeric(df["POS"], errors="coerce")
    df["SCORE"] = pd.to_numeric(df["SCORE"], errors="coerce")
    if "CRIT" in df.columns:
        df["CRIT"] = pd.to_numeric(df["CRIT"], errors="coerce")

    df["ABS_SCORE"] = df["SCORE"].abs()
    df = df.dropna(subset=["POS", "SCORE"])

    rows = []
    for chrom in df["CHROM"].dropna().unique():
        sub = df[df["CHROM"] == chrom].copy()
        if sub.empty:
            continue

        maxpos = int(sub["POS"].max())
        starts = range(1, maxpos + 1, args.step)

        for s in starts:
            e = s + args.window - 1
            win = sub[(sub["POS"] >= s) & (sub["POS"] <= e)]
            if win.empty:
                continue

            row = {
                "CHROM": chrom,
                "start": s,
                "end": e,
                "iHS_mean": win["SCORE"].mean(),
                "iHS_mean_abs": win["ABS_SCORE"].mean(),
                "iHS_max": win["SCORE"].max(),
                "iHS_min": win["SCORE"].min(),
                "iHS_max_abs": win["ABS_SCORE"].max(),
                "n_snp": win.shape[0]
            }

            if "CRIT" in win.columns:
                row["n_crit"] = int((win["CRIT"] == 1).sum())
                row["crit_frac"] = (win["CRIT"] == 1).mean()

            rows.append(row)

    out = pd.DataFrame(rows)
    out.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
