#################################################################################################
#1. Nucleotide diversity (π), the fixation index (FST), and absolute sequence divergence dXY
#################################################################################################

# 下面是run_pixy_parallel.sh的内容
#!/usr/bin/env bash
set -euo pipefail

# =========================
# 用户参数
# =========================
VCF="FMdeer.novar.noindel.filtered.vcf.gz"
POP="pop.info"
GENOME="genome.bed"

OUTDIR="/data/xb/FMdeer/04.finalSNP/FstPi"
PREFIX="FMdeer_pixy"

WINDOW=50000
STEP=10000

# 64核服务器推荐起点
JOBS=16
PIXCORES=4

# 是否保留末端不足50kb窗口：yes/no
KEEP_PARTIAL="yes"

# =========================
# 目录
# =========================
BEDDIR="${OUTDIR}/beds"
SUBBEDDIR="${BEDDIR}/by_chr"
TMPDIR="${OUTDIR}/pixy_by_chr"
MERGEDIR="${OUTDIR}/merged"
LOGDIR="${OUTDIR}/logs"

mkdir -p "${OUTDIR}" "${BEDDIR}" "${SUBBEDDIR}" "${TMPDIR}" "${MERGEDIR}" "${LOGDIR}"

ALLBED="${BEDDIR}/${PREFIX}.50kb_10kb.windows.bed"
CHR_LIST="${BEDDIR}/chrom.list"

# =========================
# 检查依赖
# =========================
for cmd in pixy parallel awk cut tail grep sort; do
    command -v "${cmd}" >/dev/null 2>&1 || {
        echo "[ERROR] Missing command: ${cmd}" >&2
        exit 1
    }
done

# =========================
# 检查输入
# =========================
for f in "${VCF}" "${POP}" "${GENOME}"; do
    [[ -s "${f}" ]] || {
        echo "[ERROR] Missing or empty file: ${f}" >&2
        exit 1
    }
done

if [[ ! -s "${VCF}.tbi" && ! -s "${VCF}.csi" ]]; then
    echo "[ERROR] Missing index for VCF: ${VCF}.tbi or ${VCF}.csi" >&2
    exit 1
fi

echo "[INFO] Input files checked."
echo "[INFO] VCF       : ${VCF}"
echo "[INFO] POP       : ${POP}"
echo "[INFO] GENOME    : ${GENOME}"
echo "[INFO] OUTDIR    : ${OUTDIR}"
echo "[INFO] WINDOW    : ${WINDOW}"
echo "[INFO] STEP      : ${STEP}"
echo "[INFO] JOBS      : ${JOBS}"
echo "[INFO] PIXCORES  : ${PIXCORES}"

# =========================
# 清理旧文件（可选）
# =========================
rm -f "${ALLBED}" "${CHR_LIST}"
rm -f "${SUBBEDDIR}"/*.bed

# =========================
# 1. 生成滑窗BED
# =========================
echo "[INFO] Generating sliding windows BED..."

if [[ "${KEEP_PARTIAL}" == "yes" ]]; then
    awk -v W="${WINDOW}" -v S="${STEP}" 'BEGIN{OFS="\t"}
    NF>=2{
        chr=$1; len=$2;
        for (start=0; start<len; start+=S) {
            end=start+W;
            if (end>len) end=len;
            print chr, start, end;
            if (end==len) break;
        }
    }' "${GENOME}" > "${ALLBED}"
else
    awk -v W="${WINDOW}" -v S="${STEP}" 'BEGIN{OFS="\t"}
    NF>=2{
        chr=$1; len=$2;
        for (start=0; start+W<=len; start+=S) {
            end=start+W;
            print chr, start, end;
        }
    }' "${GENOME}" > "${ALLBED}"
fi

echo "[INFO] BED created: ${ALLBED}"

# =========================
# 2. 染色体列表
# =========================
cut -f1 "${GENOME}" > "${CHR_LIST}"
echo "[INFO] Chromosome list written: ${CHR_LIST}"

# =========================
# 3. 按染色体拆分BED
# =========================
echo "[INFO] Splitting BED by chromosome..."

awk -v outdir="${SUBBEDDIR}" '{
    file = outdir "/" $1 ".bed";
    print $0 >> file
}' "${ALLBED}"

echo "[INFO] BED split complete."

# =========================
# 4. 并行运行pixy
# =========================
echo "[INFO] Starting parallel pixy jobs..."

parallel --joblog "${LOGDIR}/parallel_joblog.txt" \
         --results "${LOGDIR}/parallel_results" \
         -j "${JOBS}" \
         bash pixy_one_chr.sh {} "${VCF}" "${POP}" "${SUBBEDDIR}" "${TMPDIR}" "${PREFIX}" "${PIXCORES}" \
         :::: "${CHR_LIST}"

echo "[INFO] All pixy jobs completed."

# =========================
# 5. 合并结果
# =========================
merge_stat () {
    stat="$1"
    outfile="${MERGEDIR}/${PREFIX}.merged_${stat}.txt"
    first=1
    > "${outfile}"

    for f in "${TMPDIR}"/*/"${PREFIX}".*_"${stat}".txt; do
        [[ -e "$f" ]] || continue
        if [[ $first -eq 1 ]]; then
            cat "$f" >> "${outfile}"
            first=0
        else
            tail -n +2 "$f" >> "${outfile}"
        fi
    done

    if [[ -s "${outfile}" ]]; then
        echo "[INFO] Merged ${stat}: ${outfile}"
    else
        echo "[WARN] No merged output for ${stat}"
    fi
}

echo "[INFO] Merging outputs..."
merge_stat pi
merge_stat fst
merge_stat dxy

echo "[INFO] Done."

##########################################################
2. # tajimaD
##########################################################

#下面是run_tajimaD_parallel.sh的内容
#!/usr/bin/env bash
set -euo pipefail

# =========================
# 用户参数
# =========================
VCF="FMdeer.novar.noindel.filtered.vcf.gz"
POP="pop.info"
WINDOW_BED="/data/xb/FMdeer/04.finalSNP/ROH/50k.10k.window.bed"
PY_SCRIPT="tajimaD_scikit_allel_bychr.py"

DOMESTIC="domestic"
WILD="wild"

OUTDIR="/data/xb/FMdeer/04.finalSNP/TajimaD"
TMPDIR="${OUTDIR}/by_chr"
LOGDIR="${OUTDIR}/logs"

# 并发任务数
JOBS=16

# 是否强制重跑全部: yes/no
FORCE_RERUN="no"

# =========================
# 建目录
# =========================
mkdir -p "${OUTDIR}" "${TMPDIR}" "${LOGDIR}"

# =========================
# 检查依赖
# =========================
for cmd in python3 parallel bcftools awk wc tail; do
    command -v "${cmd}" >/dev/null 2>&1 || {
        echo "[ERROR] Missing command: ${cmd}" >&2
        exit 1
    }
done

# =========================
# 检查输入文件
# =========================
for f in "${VCF}" "${POP}" "${WINDOW_BED}" "${PY_SCRIPT}"; do
    [[ -s "${f}" ]] || {
        echo "[ERROR] Missing or empty file: ${f}" >&2
        exit 1
    }
done

if [[ ! -s "${VCF}.tbi" && ! -s "${VCF}.csi" ]]; then
    echo "[ERROR] Missing VCF index: ${VCF}.tbi or ${VCF}.csi" >&2
    exit 1
fi

# =========================
# 获取染色体列表
# =========================
CHR_LIST_ALL="${OUTDIR}/chrom.list.all"
CHR_LIST_TODO="${OUTDIR}/chrom.list.todo"

bcftools index -s "${VCF}" | cut -f1 > "${CHR_LIST_ALL}"

echo "[INFO] All chromosomes:"
cat "${CHR_LIST_ALL}"

# =========================
# 判断哪些染色体需要计算
# 规则：
# 1. FORCE_RERUN=yes -> 全部重跑
# 2. 输出文件不存在/为空/少于2行 -> 重跑
# =========================
> "${CHR_LIST_TODO}"

while read -r chr; do
    [[ -n "${chr}" ]] || continue
    outfile="${TMPDIR}/${chr}.tajimaD.tsv"

    if [[ "${FORCE_RERUN}" == "yes" ]]; then
        echo "${chr}" >> "${CHR_LIST_TODO}"
        continue
    fi

    if [[ ! -s "${outfile}" ]]; then
        echo "${chr}" >> "${CHR_LIST_TODO}"
        continue
    fi

    nlines=$(wc -l < "${outfile}")
    if [[ "${nlines}" -lt 2 ]]; then
        echo "${chr}" >> "${CHR_LIST_TODO}"
        continue
    fi

done < "${CHR_LIST_ALL}"

TODO_N=$(wc -l < "${CHR_LIST_TODO}" | awk '{print $1}')
ALL_N=$(wc -l < "${CHR_LIST_ALL}" | awk '{print $1}')

echo "[INFO] Total chromosomes   : ${ALL_N}"
echo "[INFO] Need to run         : ${TODO_N}"
echo "[INFO] Skip finished files : $((ALL_N - TODO_N))"

if [[ "${TODO_N}" -eq 0 ]]; then
    echo "[INFO] Nothing to run. All chromosome outputs already exist."
else
    echo "[INFO] Chromosomes to run:"
    cat "${CHR_LIST_TODO}"

    echo "[INFO] Start GNU parallel..."
    echo "[INFO] JOBS=${JOBS}"

    parallel --joblog "${LOGDIR}/parallel_joblog.txt" \
             --results "${LOGDIR}/parallel_results" \
             -j "${JOBS}" \
             python3 "${PY_SCRIPT}" \
               -v "${VCF}" \
               -p "${POP}" \
               --window-bed "${WINDOW_BED}" \
               --domestic-name "${DOMESTIC}" \
               --wild-name "${WILD}" \
               --chrom {} \
               -o "${TMPDIR}/{}.tajimaD.tsv" \
             :::: "${CHR_LIST_TODO}"

    echo "[INFO] Parallel jobs finished."
fi

# =========================
# 合并输出
# 仅合并有效文件（>=2行）
# =========================
MERGED="${OUTDIR}/tajimaD.windowed.tsv"
first=1
> "${MERGED}"

while read -r chr; do
    [[ -n "${chr}" ]] || continue
    f="${TMPDIR}/${chr}.tajimaD.tsv"

    if [[ ! -s "${f}" ]]; then
        echo "[WARN] Missing output, skip merge: ${f}" >&2
        continue
    fi

    nlines=$(wc -l < "${f}")
    if [[ "${nlines}" -lt 2 ]]; then
        echo "[WARN] Invalid output (<2 lines), skip merge: ${f}" >&2
        continue
    fi

    if [[ $first -eq 1 ]]; then
        cat "${f}" >> "${MERGED}"
        first=0
    else
        tail -n +2 "${f}" >> "${MERGED}"
    fi
done < "${CHR_LIST_ALL}"

if [[ ! -s "${MERGED}" ]]; then
    echo "[ERROR] Merged file is empty." >&2
    exit 1
fi

echo "[INFO] Merged output: ${MERGED}"

# =========================
# 检查失败/缺失染色体
# =========================
FAILED_LIST="${OUTDIR}/chrom.failed.list"
> "${FAILED_LIST}"

while read -r chr; do
    [[ -n "${chr}" ]] || continue
    f="${TMPDIR}/${chr}.tajimaD.tsv"

    if [[ ! -s "${f}" ]]; then
        echo "${chr}" >> "${FAILED_LIST}"
        continue
    fi

    nlines=$(wc -l < "${f}")
    if [[ "${nlines}" -lt 2 ]]; then
        echo "${chr}" >> "${FAILED_LIST}"
        continue
    fi
done < "${CHR_LIST_ALL}"

FAILED_N=$(wc -l < "${FAILED_LIST}" | awk '{print $1}')
if [[ "${FAILED_N}" -gt 0 ]]; then
    echo "[WARN] ${FAILED_N} chromosome(s) failed or incomplete:"
    cat "${FAILED_LIST}"
else
    echo "[INFO] All chromosome outputs look complete."
fi

echo "[INFO] Done."

