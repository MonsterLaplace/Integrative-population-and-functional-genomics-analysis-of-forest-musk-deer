#!/bin/bash
set -euo pipefail

echo "argc=$#"
echo "argv=$@"

if [ $# -ne 7 ]; then
    echo "Usage: bash 01_wgs_fastq_to_bam.sh <reference.fa> <sample1,sample2,...> <fastqdir> <bamdir> <outdir> <aligner_threads> <samtools_threads>"
    exit 1
fi

reference="$1"
samples="$2"
fastqdir="$3"
bamdir="$4"
outdir="$5"
aligner_threads="$6"
samtools_threads="$7"

logdir="${outdir}/logs/step1"
statusdir="${outdir}/status/step1"
mkdir -p "$bamdir" "$logdir" "$statusdir"

if [ ! -f "$reference" ]; then
    echo "ERROR: reference not found: $reference"
    exit 1
fi

if ! [[ "$aligner_threads" =~ ^[0-9]+$ ]]; then
    echo "ERROR: aligner_threads must be a positive integer: $aligner_threads"
    exit 1
fi

if ! [[ "$samtools_threads" =~ ^[0-9]+$ ]]; then
    echo "ERROR: samtools_threads must be a positive integer: $samtools_threads"
    exit 1
fi

samples=$(echo "$samples" | tr ',' ' ')

ref_prefix="${reference%.*}"
dict="${ref_prefix}.dict"

echo "reference=$reference"
echo "samples=$samples"
echo "fastqdir=$fastqdir"
echo "bamdir=$bamdir"
echo "outdir=$outdir"
echo "aligner_threads=$aligner_threads"
echo "samtools_threads=$samtools_threads"

echo "----------------------------------------"
echo "check/build reference index"

# minibwa index
if [ ! -f "${reference}.l2b" ] || [ ! -f "${reference}.mbw" ]; then
    echo "minibwa index not found, building..."
    time minibwa index -t "$aligner_threads" "$reference"
else
    echo "minibwa index exists"
fi

# samtools faidx
if [ ! -f "${reference}.fai" ]; then
    echo "fasta index not found, building..."
    time samtools faidx "$reference"
else
    echo "fasta index exists"
fi

# gatk sequence dictionary
if [ ! -f "$dict" ]; then
    echo "sequence dictionary not found, building..."
    time gatk CreateSequenceDictionary -R "$reference" -O "$dict"
else
    echo "sequence dictionary exists"
fi

for sample in $samples; do
    fq1="${fastqdir}/${sample}_1.fastq.gz"
    fq2="${fastqdir}/${sample}_2.fastq.gz"

    sampledir="${bamdir}/${sample}"
    qcdir="${sampledir}/qc"
    mkdir -p "$sampledir" "$qcdir"

    clean_fq1="${sampledir}/${sample}.clean.1.fastq.gz"
    clean_fq2="${sampledir}/${sample}.clean.2.fastq.gz"
    sort_bam="${sampledir}/${sample}.sort.bam"
    rg_bam="${sampledir}/${sample}.rg.bam"
    markdup_bam="${sampledir}/${sample}.markdup.bam"
    final_bam="${bamdir}/${sample}.sort.rmdup.bam"
    final_bai="${final_bam}.bai"
    metrics="${sampledir}/${sample}.markdup.metrics.txt"

    log_file="${logdir}/${sample}.log"
    done_file="${statusdir}/${sample}.done"

    echo "----------------------------------------"
    echo "sample=$sample"
    echo "fq1=$fq1"
    echo "fq2=$fq2"
    echo "final_bam=$final_bam"

    if [ -f "$final_bam" ] && [ -f "$final_bai" ] && [ -f "$done_file" ] && [ -f "$log_file" ] && grep -q "STEP1_SUCCESS" "$log_file"; then
        echo "** ${sample} step1 finished, skip **"
        continue
    fi

    rm -f "$done_file"
    rm -f "$final_bam" "$final_bai"
    rm -f "$sort_bam" "${sort_bam}.bai"
    rm -f "$rg_bam" "$markdup_bam"
    rm -f "$clean_fq1" "$clean_fq2"

    {
        echo "[$(date '+%F %T')] sample=$sample start"
        echo "fq1=$fq1"
        echo "fq2=$fq2"

        if [ ! -f "$fq1" ]; then
            echo "ERROR: fastq not found: $fq1"
            exit 1
        fi

        if [ ! -f "$fq2" ]; then
            echo "ERROR: fastq not found: $fq2"
            exit 1
        fi

        echo "raw fastqc..."
        fastqc -t "$samtools_threads" -o "$qcdir" "$fq1" "$fq2"

        echo "fastp..."
        fastp \
            -i "$fq1" \
            -I "$fq2" \
            -o "$clean_fq1" \
            -O "$clean_fq2" \
            -h "${qcdir}/${sample}.fastp.html" \
            -j "${qcdir}/${sample}.fastp.json" \
            -W 5 \
            -M 20 \
            -5 \
            -q 15 \
            -u 40 \
            -n 0 \
            -l 30 \
            -w 8

        if [ ! -s "$clean_fq1" ]; then
            echo "ERROR: fastp output file missing or empty: $clean_fq1"
            exit 1
        fi

        if [ ! -s "$clean_fq2" ]; then
            echo "ERROR: fastp output file missing or empty: $clean_fq2"
            exit 1
        fi

        if ! gzip -t "$clean_fq1"; then
            echo "ERROR: corrupted gzip file from fastp: $clean_fq1"
            exit 1
        fi

        if ! gzip -t "$clean_fq2"; then
            echo "ERROR: corrupted gzip file from fastp: $clean_fq2"
            exit 1
        fi

        echo "clean fastqc..."
        fastqc -t "$samtools_threads" -o "$qcdir" "$clean_fq1" "$clean_fq2"

        echo "multiqc..."
        multiqc "$qcdir" -o "$qcdir"

        echo "minibwa map | samtools view | samtools sort ..."
        minibwa map -t "$aligner_threads" "$reference" "$clean_fq1" "$clean_fq2" \
            | samtools view -@ "$samtools_threads" -bS - \
            | samtools sort -@ "$samtools_threads" -o "$sort_bam" -

        echo "samtools index..."
        samtools index -@ "$samtools_threads" "$sort_bam"

        echo "gatk AddOrReplaceReadGroups..."
        gatk AddOrReplaceReadGroups \
            -I "$sort_bam" \
            -O "$rg_bam" \
            -RGID "$sample" \
            -RGLB "$sample" \
            -RGPL "ILLUMINA" \
            -RGPU "${sample}.unit1" \
            -RGSM "$sample"

        echo "gatk MarkDuplicates..."
        gatk MarkDuplicates \
            -I "$rg_bam" \
            -O "$markdup_bam" \
            -M "$metrics" \
            --REMOVE_DUPLICATES true

        echo "samtools sort rmdup bam..."
        samtools sort -@ "$samtools_threads" -o "$final_bam" "$markdup_bam"

        echo "samtools index final bam..."
        samtools index -@ "$samtools_threads" "$final_bam"

        echo "remove intermediate files..."
        rm -f "$clean_fq1" "$clean_fq2"
        rm -f "$sort_bam" "${sort_bam}.bai"
        rm -f "$rg_bam"
        rm -f "$markdup_bam"

        echo "STEP1_SUCCESS"
        echo "[$(date '+%F %T')] sample=$sample done"
    } > "$log_file" 2>&1

    echo "OK $(date '+%F %T')" > "$done_file"

    echo "** ${sample}.sort.rmdup.bam done **"
done
