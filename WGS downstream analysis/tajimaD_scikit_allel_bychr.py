#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import argparse
import subprocess
import numpy as np
import pandas as pd
import allel


def parse_args():
    p = argparse.ArgumentParser(
        description="Calculate Tajima's D for two populations using predefined windows from a BED file."
    )
    p.add_argument("-v", "--vcf", required=True, help="Input bgzipped VCF/VCF.GZ")
    p.add_argument("-p", "--popfile", required=True, help="Population file: sample<TAB>group")
    p.add_argument("-o", "--out", required=True, help="Output TSV")
    p.add_argument(
        "--window-bed",
        required=True,
        help="BED file with windows: chrom start end"
    )
    p.add_argument("--domestic-name", default="domestic", help="Domestic group name")
    p.add_argument("--wild-name", default="wild", help="Wild group name")
    p.add_argument("--bcftools", default="bcftools", help="Path to bcftools executable")
    p.add_argument("--chrom", default=None, help="Only process one chromosome")
    return p.parse_args()


def check_exists(path, desc):
    if not os.path.exists(path):
        raise FileNotFoundError(f"{desc}不存在: {path}")


def run_cmd(cmd):
    r = subprocess.run(cmd, shell=True, text=True, capture_output=True)
    if r.returncode != 0:
        raise RuntimeError(f"命令执行失败:\n{cmd}\n\nstdout:\n{r.stdout}\n\nstderr:\n{r.stderr}")
    return r.stdout.strip()


def get_vcf_samples(vcf, bcftools="bcftools"):
    cmd = f"{bcftools} query -l {vcf}"
    out = run_cmd(cmd)
    if not out:
        raise ValueError("无法从VCF读取样本名")
    return np.array(out.splitlines(), dtype=object)


def get_chrom_list(vcf, bcftools="bcftools"):
    cmd = f"{bcftools} index -s {vcf}"
    out = run_cmd(cmd)
    chroms = []
    for line in out.splitlines():
        if not line.strip():
            continue
        chrom = line.split("\t")[0]
        chroms.append(chrom)
    if not chroms:
        raise ValueError("无法从VCF索引中获取染色体列表，请确认VCF已建立索引(.tbi/.csi)")
    return chroms


def load_popfile(popfile):
    pop = pd.read_csv(popfile, sep="\t", header=None, names=["sample", "group"])
    if pop.empty:
        raise ValueError("群体文件为空")
    return pop


def load_windows(window_bed):
    win = pd.read_csv(window_bed, sep="\t", header=None, names=["chrom", "start", "end"])
    if win.empty:
        raise ValueError("窗口文件为空")
    win["start"] = win["start"].astype(int)
    win["end"] = win["end"].astype(int)
    return win


def read_region_callset(vcf, region):
    callset = allel.read_vcf(
        vcf,
        region=region,
        fields=["variants/CHROM", "variants/POS", "calldata/GT"]
    )
    return callset


def calc_tajima_for_window(pos, ac, start, end):
    # BED为0-based half-open [start, end)
    # VCF POS为1-based
    # 因此窗口内条件写为: start < POS <= end
    mask = (pos > start) & (pos <= end)

    if mask.sum() < 2:
        return np.nan

    try:
        return allel.tajima_d(ac[mask], pos=pos[mask])
    except Exception:
        return np.nan


def main():
    args = parse_args()

    check_exists(args.vcf, "VCF文件")
    check_exists(args.popfile, "群体文件")
    check_exists(args.window_bed, "窗口BED文件")

    if not (os.path.exists(args.vcf + ".tbi") or os.path.exists(args.vcf + ".csi")):
        raise FileNotFoundError("VCF缺少索引文件(.tbi或.csi)，请先建立索引，例如: bcftools index -t your.vcf.gz")

    print("读取群体文件...")
    pop = load_popfile(args.popfile)

    print("读取窗口文件...")
    windows = load_windows(args.window_bed)
    print(f"窗口数: {len(windows)}")

    print("读取VCF样本名...")
    samples = get_vcf_samples(args.vcf, args.bcftools)
    vcf_sample_set = set(samples.tolist())

    domestic_samples = set(pop.loc[pop["group"] == args.domestic_name, "sample"])
    wild_samples = set(pop.loc[pop["group"] == args.wild_name, "sample"])

    if not domestic_samples:
        raise ValueError(f"群体文件中未找到群体: {args.domestic_name}")
    if not wild_samples:
        raise ValueError(f"群体文件中未找到群体: {args.wild_name}")

    domestic_idx = [i for i, s in enumerate(samples) if s in domestic_samples]
    wild_idx = [i for i, s in enumerate(samples) if s in wild_samples]

    print(f"VCF样本数: {len(samples)}")
    print(f"{args.domestic_name} 样本数(VCF中匹配到): {len(domestic_idx)}")
    print(f"{args.wild_name} 样本数(VCF中匹配到): {len(wild_idx)}")

    if len(domestic_idx) == 0:
        raise ValueError(f"VCF中没有找到 {args.domestic_name} 群体样本")
    if len(wild_idx) == 0:
        raise ValueError(f"VCF中没有找到 {args.wild_name} 群体样本")

    missing_dom = domestic_samples - vcf_sample_set
    missing_wild = wild_samples - vcf_sample_set
    if missing_dom:
        print(f"警告: {args.domestic_name} 中有 {len(missing_dom)} 个样本不在VCF中")
    if missing_wild:
        print(f"警告: {args.wild_name} 中有 {len(missing_wild)} 个样本不在VCF中")

    print("读取染色体列表...")
    chrom_list = get_chrom_list(args.vcf, args.bcftools)

    if args.chrom is not None:
        chrom_list = [args.chrom]

    print(f"待处理染色体数: {len(chrom_list)}")

    results = []

    for chrom in chrom_list:
        win_chr = windows[windows["chrom"] == chrom]
        if win_chr.empty:
            print(f"跳过 {chrom}: 窗口文件中无该染色体")
            continue

        print(f"处理: {chrom}，窗口数: {len(win_chr)}")

        try:
            callset = read_region_callset(args.vcf, chrom)
        except Exception as e:
            print(f"跳过 {chrom}: 读取失败: {e}")
            continue

        if callset is None:
            print(f"跳过 {chrom}: 无数据")
            continue

        if "variants/POS" not in callset or "calldata/GT" not in callset:
            print(f"跳过 {chrom}: 缺少POS或GT")
            continue

        pos = callset["variants/POS"]
        gt_raw = callset["calldata/GT"]

        if pos is None or gt_raw is None or len(pos) == 0:
            print(f"跳过 {chrom}: 无变异位点")
            continue

        if len(pos) < 2:
            print(f"跳过 {chrom}: 位点数不足(<2)")
            continue

        gt = allel.GenotypeArray(gt_raw)

        try:
            gt_dom = gt.take(domestic_idx, axis=1)
            gt_wild = gt.take(wild_idx, axis=1)
        except Exception as e:
            print(f"跳过 {chrom}: 提取群体样本失败: {e}")
            continue

        ac_dom = gt_dom.count_alleles()
        ac_wild = gt_wild.count_alleles()

        for _, row in win_chr.iterrows():
            start = int(row["start"])
            end = int(row["end"])

            tajd_dom = calc_tajima_for_window(pos, ac_dom, start, end)
            tajd_wild = calc_tajima_for_window(pos, ac_wild, start, end)

            results.append([
                chrom,
                start,
                end,
                tajd_dom,
                tajd_wild
            ])

    if not results:
        raise ValueError("没有生成任何结果，请检查VCF、索引、群体文件、窗口文件和群体名称")

    res = pd.DataFrame(
        results,
        columns=[
            "chrom",
            "start",
            "end",
            f"tajimaD_{args.domestic_name}",
            f"tajimaD_{args.wild_name}"
        ]
    )

    res = res.sort_values(["chrom", "start", "end"])
    res.to_csv(args.out, sep="\t", index=False)
    print(f"完成，结果已输出到: {args.out}")


if __name__ == "__main__":
    main()
