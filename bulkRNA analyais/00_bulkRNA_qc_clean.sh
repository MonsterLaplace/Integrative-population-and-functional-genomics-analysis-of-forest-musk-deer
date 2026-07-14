#!/usr/bin/env bash
set -euo pipefail
source 09.scripts/config.sh

mkdir -p 05.bulkRNA/{01.qc,02.clean}

tail -n +2 00.metadata/bulkRNA_samples.tsv | while IFS=$'\t' read -r sample r1 r2 tissue age_group; do
  echo "[INFO] bulkRNA QC ${sample}"

  ${FASTQC} -t 8 -o 05.bulkRNA/01.qc ${r1} ${r2}

  ${FASTP} \
    -i ${r1} -I ${r2} \
    -o 05.bulkRNA/02.clean/${sample}_R1.clean.fq.gz \
    -O 05.bulkRNA/02.clean/${sample}_R2.clean.fq.gz \
    -q 20 -u 30 -n 5 -l 50 -w 8 \
    -h 05.bulkRNA/02.clean/${sample}.fastp.html \
    -j 05.bulkRNA/02.clean/${sample}.fastp.json
done

${MULTIQC} 05.bulkRNA/01.qc 05.bulkRNA/02.clean -o 05.bulkRNA/01.qc/multiqc

