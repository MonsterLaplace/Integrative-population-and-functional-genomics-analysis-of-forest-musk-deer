#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd

# 输入输出文件
input_file = "all_samples.roh.hom"
output_file = "roh_count_by_sample.txt"

# 读取PLINK .hom文件
# 默认按任意空白符分隔
df = pd.read_csv(input_file, delim_whitespace=True)

# 检查列名
required_cols = ["IID", "KB"]
for col in required_cols:
    if col not in df.columns:
        raise ValueError(f"输入文件缺少必要列: {col}")

# 定义分类函数
def classify_roh(kb):
    if kb < 200:
        return "ROH_<200kb"
    elif kb <= 2500:
        return "ROH_200kb_2.5Mb"
    else:
        return "ROH_>2.5Mb"

# 添加类别列
df["ROH_Class"] = df["KB"].apply(classify_roh)

# 按样本和类别统计数量
count_table = df.groupby(["IID", "ROH_Class"]).size().unstack(fill_value=0)

# 确保三类都存在
for col in ["ROH_<200kb", "ROH_200kb_2.5Mb", "ROH_>2.5Mb"]:
    if col not in count_table.columns:
        count_table[col] = 0

# 按固定顺序排列列
count_table = count_table[["ROH_<200kb", "ROH_200kb_2.5Mb", "ROH_>2.5Mb"]]

# 添加总数列
count_table["Total_ROH"] = count_table.sum(axis=1)

# 输出
count_table = count_table.reset_index()
count_table.to_csv(output_file, sep="\t", index=False)

print(f"统计完成，结果已保存到: {output_file}")
