#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pandas as pd

CLASS_PRIORITY = {
    'LoF': 4,
    'missense_severe': 3,
    'missense_all': 2,
    'synonymous': 1
}

def pick_best_record(group):
    # 为每条记录赋优先级
    group = group.copy()
    group['priority'] = group['class'].map(CLASS_PRIORITY).fillna(0)

    # Grantham缺失时设为-1，便于排序
    if 'Grantham' in group.columns:
        group['Grantham_num'] = pd.to_numeric(group['Grantham'], errors='coerce').fillna(-1)
    else:
        group['Grantham_num'] = -1

    # 排序规则：
    # 1. class优先级高
    # 2. impact优先级：HIGH > MODERATE > LOW > MODIFIER
    impact_priority = {'HIGH':4, 'MODERATE':3, 'LOW':2, 'MODIFIER':1}
    if 'impact' in group.columns:
        group['impact_priority'] = group['impact'].map(impact_priority).fillna(0)
    else:
        group['impact_priority'] = 0

    group = group.sort_values(
        by=['priority', 'impact_priority', 'Grantham_num'],
        ascending=[False, False, False]
    )

    best = group.iloc[0].copy()

    # 附加信息：记录该位点原始有多少条转录本注释
    best['n_annotations'] = group.shape[0]

    # 记录该位点涉及到哪些类别
    best['all_classes'] = ",".join(sorted(group['class'].dropna().unique().tolist()))

    # 记录涉及哪些effect
    if 'effect' in group.columns:
        best['all_effects'] = ",".join(sorted(group['effect'].dropna().unique().tolist()))
    else:
        best['all_effects'] = ""

    # 记录涉及哪些gene
    if 'gene' in group.columns:
        genes = [g for g in group['gene'].dropna().astype(str).unique().tolist() if g not in ['', '.']]
        best['all_genes'] = ",".join(sorted(genes))
    else:
        best['all_genes'] = ""

    # 记录涉及哪些transcript
    if 'transcript' in group.columns:
        txs = [t for t in group['transcript'].dropna().astype(str).unique().tolist() if t not in ['', '.']]
        best['all_transcripts'] = ",".join(sorted(txs))
    else:
        best['all_transcripts'] = ""

    # 去掉中间列
    best = best.drop(labels=['priority', 'impact_priority', 'Grantham_num'], errors='ignore')
    return best

def main():
    parser = argparse.ArgumentParser(description="Collapse transcript-level functional variants to site-level.")
    parser.add_argument("-i", "--input", required=True, help="Input transcript-level TSV (functional_variants_with_grantham.tsv)")
    parser.add_argument("-o", "--output", required=True, help="Output site-level TSV")
    args = parser.parse_args()

    df = pd.read_csv(args.input, sep="\t", dtype=str)

    required_cols = ['CHROM', 'POS', 'REF', 'ALT', 'class']
    for col in required_cols:
        if col not in df.columns:
            raise ValueError(f"Required column missing: {col}")

    grouped = df.groupby(['CHROM', 'POS', 'REF', 'ALT'], as_index=False, group_keys=False)
    out = grouped.apply(pick_best_record)

    # 排序
    try:
        out['POS_num'] = pd.to_numeric(out['POS'], errors='coerce')
        out = out.sort_values(by=['CHROM', 'POS_num', 'REF', 'ALT']).drop(columns=['POS_num'])
    except:
        out = out.sort_values(by=['CHROM', 'POS', 'REF', 'ALT'])

    out.to_csv(args.output, sep="\t", index=False)

    # summary
    print("=== Site-level summary ===")
    print(out['class'].value_counts(dropna=False).to_string())
    print(f"\nOutput written to: {args.output}")

if __name__ == "__main__":
    main()
