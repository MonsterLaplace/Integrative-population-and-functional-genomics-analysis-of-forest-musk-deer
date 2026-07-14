#!/usr/bin/env bash
set -euo pipefail
source 09.scripts/config.sh

mkdir -p 05.bulkRNA/{03.bam,04.counts}

tail -n +2 00.metadata/bulkRNA_samples.tsv | while IFS=$'\t' read -r sample r1 r2 tissue age_group; do
  bam=05.bulkRNA/03.bam/${sample}.bam
  bai=${bam}.bai

  if [[ -s "${bam}" && -s "${bai}" ]]; then
    echo "[INFO] Skip ${sample}: ${bam} and ${bai} already exist"
    continue
  fi

  if [[ -s "${bam}" && ! -s "${bai}" ]]; then
    echo "[INFO] Index existing BAM for ${sample}"
    ${SAMTOOLS} index "${bam}"
    continue
  fi

  echo "[INFO] STAR ${sample}"

  ${STAR} --runThreadN 64 \
    --genomeDir ${PROJECT_DIR}/01.reference/star_index \
    --readFilesIn 05.bulkRNA/02.clean/${sample}_R1.clean.fq.gz 05.bulkRNA/02.clean/${sample}_R2.clean.fq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix 05.bulkRNA/03.bam/${sample}.

  mv 05.bulkRNA/03.bam/${sample}.Aligned.sortedByCoord.out.bam "${bam}"
  ${SAMTOOLS} index "${bam}"
done

${FEATURECOUNTS} -T 16 -p \
  -a ${REF_GTF} \
  -t exon -g gene_id \
  -o 05.bulkRNA/04.counts/gene_counts.txt \
  05.bulkRNA/03.bam/*.bam

