#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
from collections import defaultdict

#=========================
# 1. 输入文件
#=========================
roh_file = "all_samples.roh.hom"
window_file = "50k.10k.window.bed"
output_file = "roh_window_frequency.tsv"

# 样本总数
sample_list = [
    "BY01","BY02","BY03","BY04","BY05","BY06","BY07","BY08","BY09",
    "EQL1","EQL2","LZ01","LZ02","LZ03",
    "MF101","MF105","MF12","MF134","MF138","MF140","MF154","MF159",
    "MF161","MF165","MF17","MF33","MF43","MF45","MF46","MF87","MF99",
    "WQL3","WQL5","wild01","wild02","wild03","wild04","wild05"
]
n_samples = len(sample_list)

#=========================
# 2. 读取文件
#=========================
roh = pd.read_csv(roh_file, delim_whitespace=True)
windows = pd.read_csv(window_file, sep="\t", header=None, names=["CHROM", "start", "end"])

# 检查窗口文件是否真的是50kb窗口
if len(windows) <= 50:
    print("警告：你的窗口文件行数非常少，当前看起来不像真正的50kb窗口文件。")
    print("请确认 window 文件是否已经按全基因组划分为连续50kb区间。")

#=========================
# 3. 染色体名称统一函数
#=========================
def convert_chr(x):
    x = str(x)
    if x.isdigit():
        return f"Chr{int(x):02d}"
    elif x.upper() == "X":
        return "ChrX"
    elif x.upper() == "Y":
        return "ChrY"
    else:
        return x

roh["CHR"] = roh["CHR"].apply(convert_chr)
windows["CHROM"] = windows["CHROM"].apply(convert_chr)

#=========================
# 4. 按染色体整理窗口
#=========================
windows_by_chr = {}
for chr_name, subdf in windows.groupby("CHROM"):
    subdf = subdf.sort_values("start").reset_index(drop=True)
    windows_by_chr[chr_name] = subdf

# 每个窗口被哪些样本覆盖
window_samples = defaultdict(set)

# 每个窗口累计重叠的ROH条数
window_roh_counts = defaultdict(int)

#=========================
# 5. 遍历ROH，寻找重叠窗口
#=========================
for _, row in roh.iterrows():
    chr_name = row["CHR"]
    sample = row["IID"]
    start = int(row["POS1"]) - 1   # 转为BED风格
    end = int(row["POS2"])

    if chr_name not in windows_by_chr:
        continue

    chr_windows = windows_by_chr[chr_name]

    # 找所有重叠窗口
    overlaps = chr_windows[(chr_windows["end"] > start) & (chr_windows["start"] < end)]

    for _, w in overlaps.iterrows():
        key = (w["CHROM"], int(w["start"]), int(w["end"]))
        window_samples[key].add(sample)
        window_roh_counts[key] += 1

#=========================
# 6. 输出结果
#=========================
result = []
for _, w in windows.iterrows():
    key = (w["CHROM"], int(w["start"]), int(w["end"]))

    covered_samples = len(window_samples[key])      # 被多少个样本覆盖
    total_roh_count = window_roh_counts[key]        # 总共重叠了多少条ROH

    roh_freq = covered_samples / n_samples
    mean_roh_count = total_roh_count / n_samples

    result.append([
        w["CHROM"],
        int(w["start"]),
        int(w["end"]),
        roh_freq,
        mean_roh_count
    ])

result_df = pd.DataFrame(result, columns=[
    "CHROM", "start", "end", "roh_freq", "mean_roh_count"
])

result_df.to_csv(output_file, sep="\t", index=False)

print(f"完成！结果已保存到: {output_file}")
