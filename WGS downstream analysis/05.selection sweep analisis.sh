##################################
#1. Fst/pi
##################################
Rscript FstvsPi.R

##################################
#2. XP-CLR
##################################
#!/bin/bash
set -euo pipefail

vcf=/data/xb/FMdeer/04.finalSNP/all_samples.phased.vcf.gz
samplesA=domestic.list
samplesB=wild.list
outdir=/data/xb/FMdeer/04.finalSNP/XPSweep/XPCLR
mkdir -p ${outdir}/logs

mkdir -p $outdir
mkdir -p $outdir/logs

for chr in Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 Chr11 Chr12 Chr13 Chr14 Chr15 Chr16 Chr17 Chr18 Chr19 Chr20 Chr21 Chr22 Chr23 Chr24 Chr25 Chr26 Chr27 Chr28 Chr29 Chr30 Chr31 ChrX ChrY
do
  echo "[$(date '+%F %T')] Start XP-CLR on $chr"

  xpclrs \
    --input $vcf \
    --out $outdir/${chr}.xpclr.tsv \
    --samplesA $samplesA \
    --samplesB $samplesB \
    --chr $chr \
    --rrate 1e-8 \
    --ld 0.95 \
    --maxsnps 200 \
    --minsnps 10 \
    --size 50000 \
    --step 10000 \
    --threads 24 \
    --format tsv \
    --log info \
    > $outdir/logs/${chr}.log 2>&1

  echo "[$(date '+%F %T')] Finished XP-CLR on $chr"
done

Rscript xpclrs_postprocess.R

##################################
#3. XP-EHH
##################################
#!/bin/bash
set -euo pipefail

# User settings
QUERY_VCF="domestic.vcf.gz"
REF_VCF="wild.vcf.gz"
THREADS=16
OUTPREFIX="domestic_vs_wild"
DO_NORM="yes"   # yes / no

# Prepare directories

mkdir -p 00.logs
mkdir -p 01.by_chr_vcf
mkdir -p 02.by_chr_map
mkdir -p 03.by_chr_pos
mkdir -p 04.by_chr_out
mkdir -p 05.merged

LOG="00.logs/xpehh_pipeline.log"
echo "XP-EHH pipeline started at $(date)" > "${LOG}"


# Check input files

for f in "${QUERY_VCF}" "${REF_VCF}"; do
    if [ ! -f "$f" ]; then
        echo "[ERROR] File not found: $f" | tee -a "${LOG}"
        exit 1
    fi
done

if [ ! -f "${QUERY_VCF}.tbi" ] && [ ! -f "${QUERY_VCF}.csi" ]; then
    echo "[INFO] Indexing ${QUERY_VCF}" | tee -a "${LOG}"
    tabix -p vcf "${QUERY_VCF}"
fi

if [ ! -f "${REF_VCF}.tbi" ] && [ ! -f "${REF_VCF}.csi" ]; then
    echo "[INFO] Indexing ${REF_VCF}" | tee -a "${LOG}"
    tabix -p vcf "${REF_VCF}"
fi


# Chromosome list

bcftools query -f '%CHROM\n' "${QUERY_VCF}" | uniq > chr.list
echo "[INFO] Chromosome list written to chr.list" | tee -a "${LOG}"


# Run selscan per chromosome

while read chr; do
    [ -z "${chr}" ] && continue
    echo "==== Processing ${chr} ====" | tee -a "${LOG}"

    QCHR="01.by_chr_vcf/query.${chr}.vcf.gz"
    RCHR="01.by_chr_vcf/ref.${chr}.vcf.gz"
    MAP="02.by_chr_map/${chr}.map"
    QPOS="03.by_chr_pos/query.${chr}.pos"
    RPOS="03.by_chr_pos/ref.${chr}.pos"
    OUT="04.by_chr_out/${OUTPREFIX}.${chr}"

    # 1. Split by chromosome
    bcftools view -r "${chr}" "${QUERY_VCF}" -Oz -o "${QCHR}"
    bcftools view -r "${chr}" "${REF_VCF}"   -Oz -o "${RCHR}"

    tabix -f -p vcf "${QCHR}"
    tabix -f -p vcf "${RCHR}"

    # 2. Skip empty chromosome
    qn=$(bcftools view -H "${QCHR}" | wc -l)
    rn=$(bcftools view -H "${RCHR}" | wc -l)

    if [ "${qn}" -eq 0 ] || [ "${rn}" -eq 0 ]; then
        echo "[WARN] ${chr}: no variants in one of the VCFs, skip." | tee -a "${LOG}"
        continue
    fi

    # 3. Compare positions
    bcftools query -f '%CHROM\t%POS\n' "${QCHR}" > "${QPOS}"
    bcftools query -f '%CHROM\t%POS\n' "${RCHR}" > "${RPOS}"

    if ! diff -q "${QPOS}" "${RPOS}" > /dev/null; then
        echo "[WARN] ${chr}: variant positions differ between query/ref, skip this chromosome." | tee -a "${LOG}"
        echo "[WARN] First differences:" | tee -a "${LOG}"
        diff "${QPOS}" "${RPOS}" | head -20 | tee -a "${LOG}"
        continue
    fi

    # 4. Build map
    bcftools query -f '%CHROM\t%POS\n' "${QCHR}" | \
    awk 'BEGIN{OFS="\t"} {print $1,$1"_"$2,$2/1000000,$2}' > "${MAP}"

    # 5. Run selscan
    echo "[INFO] Running selscan for ${chr}" | tee -a "${LOG}"
    selscan --xpehh \
      --vcf "${QCHR}" \
      --vcf-ref "${RCHR}" \
      --map "${MAP}" \
      --out "${OUT}" \
      --threads "${THREADS}" \
      2>&1 | tee -a "${LOG}"

done < chr.list

# Merge raw outputs

echo "[INFO] Merging raw XP-EHH outputs..." | tee -a "${LOG}"

mapfile -t RAW_FILES < <(find 04.by_chr_out -type f -name "*.xpehh.out" | sort)

if [ "${#RAW_FILES[@]}" -gt 0 ]; then
    cat "${RAW_FILES[0]}" > "05.merged/${OUTPREFIX}.genome.xpehh.out"
    for f in "${RAW_FILES[@]:1}"; do
        tail -n +2 "${f}" >> "05.merged/${OUTPREFIX}.genome.xpehh.out"
    done
    echo "[INFO] Merged raw file: 05.merged/${OUTPREFIX}.genome.xpehh.out" | tee -a "${LOG}"
else
    echo "[WARN] No raw XP-EHH output files found." | tee -a "${LOG}"
fi


# Joint normalization using selscan norm

if [ "${DO_NORM}" = "yes" ]; then
    if selscan norm --help >/dev/null 2>&1; then
        echo "[INFO] Running joint normalization with selscan norm..." | tee -a "${LOG}"

        if [ "${#RAW_FILES[@]}" -eq 0 ]; then
            echo "[WARN] No raw files available for normalization, skip norm." | tee -a "${LOG}"
        else
            selscan norm --xpehh --files "${RAW_FILES[@]}" 2>&1 | tee -a "${LOG}"

            mapfile -t NORM_FILES < <(find 04.by_chr_out -type f -name "*.xpehh.out.norm" | sort)

            if [ "${#NORM_FILES[@]}" -gt 0 ]; then
                cat "${NORM_FILES[0]}" > "05.merged/${OUTPREFIX}.genome.xpehh.out.norm"
                for f in "${NORM_FILES[@]:1}"; do
                    tail -n +2 "${f}" >> "05.merged/${OUTPREFIX}.genome.xpehh.out.norm"
                done
                echo "[INFO] Merged normed file: 05.merged/${OUTPREFIX}.genome.xpehh.out.norm" | tee -a "${LOG}"
            else
                echo "[WARN] No .norm files found after normalization." | tee -a "${LOG}"
            fi
        fi
    else
        echo "[WARN] selscan norm not available, skip normalization." | tee -a "${LOG}"
    fi
fi

echo "XP-EHH pipeline finished at $(date)" | tee -a "${LOG}"

##################################
#34. iHS
##################################
#!/bin/bash
set -euo pipefail

# User settings

VCF="domestic.vcf.gz"
THREADS=16
OUTPREFIX="domestic"
DO_NORM="yes"   # yes / no


# Prepare directories

mkdir -p 00.logs
mkdir -p 11.ihs_by_chr_vcf
mkdir -p 12.ihs_by_chr_map
mkdir -p 13.ihs_by_chr_out
mkdir -p 14.ihs_merged

LOG="00.logs/ihs_pipeline.log"
echo "iHS pipeline started at $(date)" > "${LOG}"


# Check input

if [ ! -f "${VCF}" ]; then
    echo "[ERROR] File not found: ${VCF}" | tee -a "${LOG}"
    exit 1
fi

if [ ! -f "${VCF}.tbi" ] && [ ! -f "${VCF}.csi" ]; then
    echo "[INFO] Indexing ${VCF}" | tee -a "${LOG}"
    tabix -p vcf "${VCF}"
fi


# Chromosome list

bcftools query -f '%CHROM\n' "${VCF}" | uniq > chr.list
echo "[INFO] Chromosome list written to chr.list" | tee -a "${LOG}"


# Run selscan per chromosome

while read chr; do
    [ -z "${chr}" ] && continue
    echo "==== Processing ${chr} ====" | tee -a "${LOG}"

    CHRVCF="11.ihs_by_chr_vcf/${OUTPREFIX}.${chr}.vcf.gz"
    MAP="12.ihs_by_chr_map/${chr}.map"
    OUT="13.ihs_by_chr_out/${OUTPREFIX}.${chr}"

    # 1. Split by chromosome
    bcftools view -r "${chr}" "${VCF}" -Oz -o "${CHRVCF}"
    tabix -f -p vcf "${CHRVCF}"

    # 2. Skip empty chromosome
    n=$(bcftools view -H "${CHRVCF}" | wc -l)
    if [ "${n}" -eq 0 ]; then
        echo "[WARN] ${chr}: no variants, skip." | tee -a "${LOG}"
        continue
    fi

    # 3. Build map
    bcftools query -f '%CHROM\t%POS\n' "${CHRVCF}" | \
    awk 'BEGIN{OFS="\t"} {print $1,$1"_"$2,$2/1000000,$2}' > "${MAP}"

    # 4. Run selscan iHS
    echo "[INFO] Running selscan --ihs for ${chr}" | tee -a "${LOG}"
    selscan --ihs \
      --vcf "${CHRVCF}" \
      --map "${MAP}" \
      --out "${OUT}" \
      --threads "${THREADS}" \
      2>&1 | tee -a "${LOG}"

done < chr.list

# Merge raw outputs

echo "[INFO] Merging raw iHS outputs..." | tee -a "${LOG}"

mapfile -t RAW_FILES < <(find 13.ihs_by_chr_out -type f -name "*.ihs.out" | sort)

if [ "${#RAW_FILES[@]}" -gt 0 ]; then
    cat "${RAW_FILES[0]}" > "14.ihs_merged/${OUTPREFIX}.genome.ihs.out"
    for f in "${RAW_FILES[@]:1}"; do
        cat "${f}" >> "14.ihs_merged/${OUTPREFIX}.genome.ihs.out"
    done
    echo "[INFO] Merged raw file: 14.ihs_merged/${OUTPREFIX}.genome.ihs.out" | tee -a "${LOG}"
else
    echo "[WARN] No raw iHS output files found." | tee -a "${LOG}"
fi


# Joint normalization using selscan norm

if [ "${DO_NORM}" = "yes" ]; then
    if selscan norm --help >/dev/null 2>&1; then
        echo "[INFO] Running joint normalization with selscan norm..." | tee -a "${LOG}"

        if [ "${#RAW_FILES[@]}" -eq 0 ]; then
            echo "[WARN] No raw files available for normalization, skip norm." | tee -a "${LOG}"
        else
            selscan norm --ihs --files "${RAW_FILES[@]}" 2>&1 | tee -a "${LOG}"

            mapfile -t NORM_FILES < <(find 13.ihs_by_chr_out -type f -name "*.ihs.out.*.norm" | sort)

            if [ "${#NORM_FILES[@]}" -gt 0 ]; then
                cat "${NORM_FILES[0]}" > "14.ihs_merged/${OUTPREFIX}.genome.ihs.out.norm"
                for f in "${NORM_FILES[@]:1}"; do
                    cat "${f}" >> "14.ihs_merged/${OUTPREFIX}.genome.ihs.out.norm"
                done
                echo "[INFO] Merged normed file: 14.ihs_merged/${OUTPREFIX}.genome.ihs.out.norm" | tee -a "${LOG}"
            else
                echo "[WARN] No .norm files found after normalization." | tee -a "${LOG}"
            fi
        fi
    else
        echo "[WARN] selscan norm not available, skip normalization." | tee -a "${LOG}"
    fi
fi

echo "iHS pipeline finished at $(date)" | tee -a "${LOG}"



