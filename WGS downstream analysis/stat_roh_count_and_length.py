#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd

input_file = "all_samples.roh.hom"
output_file = "roh_count_length_by_sample.txt"

sample_order = [
    "BY01","BY02","BY03","BY04","BY05","BY06","BY07","BY08","BY09",
    "EQL1","EQL2",
    "LZ01","LZ02","LZ03",
    "MF101","MF105","MF12","MF134","MF138","MF140","MF154","MF159","MF161","MF165",
    "MF17","MF33","MF43","MF45","MF46","MF87","MF99",
    "WQL3","WQL5",
    "wild01","wild02","wild03","wild04","wild05"
]

df = pd.read_csv(input_file, delim_whitespace=True)

required_cols = ["IID", "KB"]
for col in required_cols:
    if col not in df.columns:
        raise ValueError(f"输入文件缺少必要列: {col}")

def classify_roh(kb):
    if kb < 200:
        return "<200kb"
    elif kb <= 2500:
        return "200kb_2.5Mb"
    else:
        return ">2.5Mb"

df["ROH_Class"] = df["KB"].apply(classify_roh)

# 统计数量
count_table = df.groupby(["IID", "ROH_Class"]).size().unstack(fill_value=0)

# 统计长度和
length_table = df.groupby(["IID", "ROH_Class"])["KB"].sum().unstack(fill_value=0)

# 补齐列
classes = ["<200kb", "200kb_2.5Mb", ">2.5Mb"]
for c in classes:
    if c not in count_table.columns:
        count_table[c] = 0
    if c not in length_table.columns:
        length_table[c] = 0

count_table = count_table[classes]
length_table = length_table[classes]

# 合并输出
result = pd.DataFrame(index=count_table.index)
result["N_<200kb"] = count_table["<200kb"]
result["KB_<200kb"] = length_table["<200kb"]
result["N_200kb_2.5Mb"] = count_table["200kb_2.5Mb"]
result["KB_200kb_2.5Mb"] = length_table["200kb_2.5Mb"]
result["N_>2.5Mb"] = count_table[">2.5Mb"]
result["KB_>2.5Mb"] = length_table[">2.5Mb"]
result["Total_N"] = result["N_<200kb"] + result["N_200kb_2.5Mb"] + result["N_>2.5Mb"]
result["Total_KB"] = result["KB_<200kb"] + result["KB_200kb_2.5Mb"] + result["KB_>2.5Mb"]

# 按样本顺序输出
result = result.reindex(sample_order, fill_value=0)
result = result.reset_index().rename(columns={"index": "IID"})

result.to_csv(output_file, sep="\t", index=False)

print(f"统计完成，结果已保存到: {output_file}")
