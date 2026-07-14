#!/bin/bash
set -euo pipefail

show_help() {
    cat << EOF
Usage:
  bash 00_run_all_wgs_pipeline.sh <datadir> <outdir> <reference.fa> <outname> <aligner_threads> <samtools_threads> <steps>

Description:
  Run the whole WGS pipeline from FASTQ to final VCF, with:
    - resumable execution
    - step selection
    - separate raw-data and result directories
    - status/log-based completion checking

Arguments:
  <datadir>
      Raw data directory.
      This directory should contain:
        1) reference genome fasta
        2) paired-end FASTQ files

      Example:
        /home/xb/raw_data

  <outdir>
      Output/result directory.
      The pipeline will create subdirectories here:
        bam/
        gvcf/
        vcf/
        logs/
        status/

      Example:
        /home/xb/result

  <reference.fa>
      Reference genome fasta file name inside <datadir>.

      Example:
        Calypte_anna.fna

      Full path resolved by script:
        <datadir>/<reference.fa>

  <outname>
      Prefix for final merged output files.

      Example:
        all_samples

      Final files may look like:
        <outdir>/vcf/all_samples.HC.snp.pass.vcf.gz
        <outdir>/vcf/all_samples.novar.noindel.filtered.vcf.gz

  <aligner_threads>
      Thread number for minibwa index/map.

      Example:
        16

  <samtools_threads>
      Thread number for samtools and some QC-related tools.

      Example:
        4

  <steps>
      Which step(s) to run.

      Allowed values:
        1         run FASTQ -> BAM only
        2         run BAM -> GVCF only
        3         run GVCF -> VCF only
        1,2       run step1 and step2
        2,3       run step2 and step3
        1,2,3     run all three steps
        all       same as 1,2,3

Pipeline Steps:
  Step 1:
      FASTQ -> BAM
      Includes:
        - minibwa index
        - samtools faidx
        - gatk CreateSequenceDictionary
        - fastqc
        - fastp
        - multiqc
        - minibwa map
        - samtools sort/index
        - gatk AddOrReplaceReadGroups
        - gatk MarkDuplicates
      Final output:
        <outdir>/bam/<sample>.sort.rmdup.bam

  Step 2:
      BAM -> GVCF
      Includes:
        - gatk HaplotypeCaller
      Final output:
        <outdir>/gvcf/<sample>.g.vcf.gz

  Step 3:
      GVCF -> final VCF
      Includes:
        - gatk CombineGVCFs
        - gatk GenotypeGVCFs
        - gatk SelectVariants
        - gatk VariantFiltration
        - bcftools processing
      Final output:
        <outdir>/vcf/<outname>.novar.noindel.filtered.vcf.gz

Parallel interval settings for Step 2/3:
  WGS_INTERVAL_PARTS
      Target number of genome chunks generated from <reference.fa>.fai
      when WGS_INTERVAL_MODE=balanced or auto selects balanced mode.
      Recommended: 10-20. Default: 20.

  WGS_INTERVAL_MODE
      Genome interval splitting strategy. Default: auto.
      auto:
        If chromosome-like contigs explain most genome length, split by chromosome
        and merge unplaced scaffolds into one final interval file.
        Otherwise split by cumulative bp size into balanced chunks.
      balanced:
        Split by cumulative bp size into similarly sized chunks.
      chromosome:
        Split chromosome-like contigs one per interval file and merge unplaced
        scaffolds into one final interval file.

  WGS_PARALLEL_JOBS
      Maximum number of interval jobs running at the same time.
      Default: 20.

  WGS_HC_THREADS
      Native PairHMM threads per HaplotypeCaller interval job.
      Used in Step 2 only. Default: 4.

  WGS_GATK_JAVA_OPTIONS
      Java options passed to GATK.
      Default: -Xmx24g.

  WGS_BCFTOOLS_THREADS
      Compression/filtering threads for bcftools commands in Step 3.
      Default: 8.

      Example:
        WGS_INTERVAL_MODE=auto WGS_INTERVAL_PARTS=20 WGS_PARALLEL_JOBS=20 WGS_HC_THREADS=4 WGS_GATK_JAVA_OPTIONS="-Xmx24g" bash 00_run_all_wgs_pipeline.sh ...

FASTQ Naming Requirement:
  Current script detects samples from:
    *_1.fastq.gz

  And expects paired files:
    <sample>_1.fastq.gz
    <sample>_2.fastq.gz

  Example:
    SRR949793_1.fastq.gz
    SRR949793_2.fastq.gz

Resume Logic:
  A step/sample is considered completed only if:
    1) result file exists
    2) index file exists (if required)
    3) done-status file exists
    4) step log contains success tag

  Success tags:
    STEP1_SUCCESS
    STEP2_SUCCESS
    STEP3_SUCCESS

Examples:
  1. Run all steps:
     bash 00_run_all_wgs_pipeline.sh /home/xb/raw_data /home/xb/result Calypte_anna.fna all_samples 16 4 all

  2. Run step1 only:
     bash 00_run_all_wgs_pipeline.sh /home/xb/raw_data /home/xb/result Calypte_anna.fna all_samples 16 4 1

  3. Run step2 and step3:
     bash 00_run_all_wgs_pipeline.sh /home/xb/raw_data /home/xb/result Calypte_anna.fna all_samples 16 4 2,3

  4. Run step3 only:
     bash 00_run_all_wgs_pipeline.sh /home/xb/raw_data /home/xb/result Calypte_anna.fna all_samples 16 4 3
EOF
}

if [ $# -eq 1 ]; then
    if [ "$1" = "-h" ] || [ "$1" = "--help" ]; then
        show_help
        exit 0
    fi
fi

echo "argc=$#"
echo "argv=$@"

if [ $# -ne 7 ]; then
    echo "ERROR: wrong number of arguments."
    echo
    show_help
    exit 1
fi

datadir="$1"
outdir="$2"
reference_name="$3"
outname="$4"
aligner_threads="$5"
samtools_threads="$6"
steps="$7"

script_dir="$(cd "$(dirname "$0")" && pwd)"

script1="${script_dir}/01_wgs_fastq_to_bam.sh"
script2="${script_dir}/02_wgs_bam_to_gvcf.sh"
script3="${script_dir}/03_wgs_gvcf_to_vcf.sh"

reference="${datadir}/${reference_name}"
fastqdir="$datadir"
bamdir="${outdir}/bam"
gvcfdir="${outdir}/gvcf"
vcfdir="${outdir}/vcf"
logdir="${outdir}/logs"
statusdir="${outdir}/status"

mkdir -p "$bamdir" "$gvcfdir" "$vcfdir"
mkdir -p "$logdir/step1" "$logdir/step2" "$logdir/step3" "$logdir/pipeline"
mkdir -p "$statusdir/step1" "$statusdir/step2" "$statusdir/step3"

pipeline_log="${logdir}/pipeline/run_all.log"

if [ ! -d "$datadir" ]; then
    echo "ERROR: datadir not found: $datadir" | tee -a "$pipeline_log"
    exit 1
fi

if [ ! -f "$reference" ]; then
    echo "ERROR: reference not found: $reference" | tee -a "$pipeline_log"
    exit 1
fi

if [ ! -f "$script1" ] || [ ! -f "$script2" ] || [ ! -f "$script3" ]; then
    echo "ERROR: one or more scripts not found in $script_dir" | tee -a "$pipeline_log"
    exit 1
fi

if ! [[ "$aligner_threads" =~ ^[0-9]+$ ]]; then
    echo "ERROR: aligner_threads must be a positive integer: $aligner_threads" | tee -a "$pipeline_log"
    exit 1
fi

if ! [[ "$samtools_threads" =~ ^[0-9]+$ ]]; then
    echo "ERROR: samtools_threads must be a positive integer: $samtools_threads" | tee -a "$pipeline_log"
    exit 1
fi

echo "script_dir=$script_dir" | tee -a "$pipeline_log"
echo "datadir=$datadir" | tee -a "$pipeline_log"
echo "outdir=$outdir" | tee -a "$pipeline_log"
echo "reference=$reference" | tee -a "$pipeline_log"
echo "fastqdir=$fastqdir" | tee -a "$pipeline_log"
echo "bamdir=$bamdir" | tee -a "$pipeline_log"
echo "gvcfdir=$gvcfdir" | tee -a "$pipeline_log"
echo "vcfdir=$vcfdir" | tee -a "$pipeline_log"
echo "outname=$outname" | tee -a "$pipeline_log"
echo "aligner_threads=$aligner_threads" | tee -a "$pipeline_log"
echo "samtools_threads=$samtools_threads" | tee -a "$pipeline_log"
echo "steps=$steps" | tee -a "$pipeline_log"

sample_list=$(find "$fastqdir" -maxdepth 1 -type f -name "*_1.fastq.gz" \
    | sed 's#.*/##' \
    | sed 's/_1.fastq.gz$//' \
    | sort)

if [ -z "$sample_list" ]; then
    echo "ERROR: no fastq files found matching *_1.fastq.gz in $fastqdir" | tee -a "$pipeline_log"
    exit 1
fi

samples=$(echo "$sample_list" | paste -sd "," -)
samples_space=$(echo "$samples" | tr ',' ' ')

echo "samples=$samples" | tee -a "$pipeline_log"

run_step1=false
run_step2=false
run_step3=false

if [ "$steps" = "all" ]; then
    run_step1=true
    run_step2=true
    run_step3=true
else
    IFS=',' read -r -a step_array <<< "$steps"
    for s in "${step_array[@]}"; do
        case "$s" in
            1) run_step1=true ;;
            2) run_step2=true ;;
            3) run_step3=true ;;
            *)
                echo "ERROR: invalid step value: $s" | tee -a "$pipeline_log"
                exit 1
                ;;
        esac
    done
fi

check_step1_sample_done() {
    local sample="$1"
    local bam="${bamdir}/${sample}.sort.rmdup.bam"
    local bai="${bam}.bai"
    local done_file="${statusdir}/step1/${sample}.done"
    local log_file="${logdir}/step1/${sample}.log"
    [ -f "$bam" ] && [ -f "$bai" ] && [ -f "$done_file" ] && [ -f "$log_file" ] && grep -q "STEP1_SUCCESS" "$log_file"
}

check_step2_sample_done() {
    local sample="$1"
    local gvcf="${gvcfdir}/${sample}.g.vcf.gz"
    local tbi="${gvcf}.tbi"
    local done_file="${statusdir}/step2/${sample}.done"
    local log_file="${logdir}/step2/${sample}.log"
    [ -f "$gvcf" ] && [ -f "$tbi" ] && [ -f "$done_file" ] && [ -f "$log_file" ] && grep -q "STEP2_SUCCESS" "$log_file"
}

check_step3_done() {
    local final_vcf="${vcfdir}/${outname}.novar.noindel.filtered.vcf.gz"
    local final_tbi="${final_vcf}.tbi"
    local done_file="${statusdir}/step3/${outname}.done"
    local log_file="${logdir}/step3/${outname}.log"
    [ -f "$final_vcf" ] && [ -f "$final_tbi" ] && [ -f "$done_file" ] && [ -f "$log_file" ] && grep -q "STEP3_SUCCESS" "$log_file"
}

step1_done=true
for sample in $samples_space; do
    if ! check_step1_sample_done "$sample"; then
        step1_done=false
        break
    fi
done

step2_done=true
for sample in $samples_space; do
    if ! check_step2_sample_done "$sample"; then
        step2_done=false
        break
    fi
done

step3_done=false
if check_step3_done; then
    step3_done=true
fi

echo "run_step1=$run_step1" | tee -a "$pipeline_log"
echo "run_step2=$run_step2" | tee -a "$pipeline_log"
echo "run_step3=$run_step3" | tee -a "$pipeline_log"
echo "step1_done=$step1_done" | tee -a "$pipeline_log"
echo "step2_done=$step2_done" | tee -a "$pipeline_log"
echo "step3_done=$step3_done" | tee -a "$pipeline_log"

if [ "$run_step1" = true ]; then
    if [ "$step1_done" = false ]; then
        echo "Step 1: FASTQ -> BAM" | tee -a "$pipeline_log"
        time bash "$script1" \
            "$reference" \
            "$samples" \
            "$fastqdir" \
            "$bamdir" \
            "$outdir" \
            "$aligner_threads" \
            "$samtools_threads"
    else
        echo "Step 1 already finished, skip." | tee -a "$pipeline_log"
    fi
fi

if [ "$run_step2" = true ]; then
    if [ "$step2_done" = false ]; then
        echo "Step 2: BAM -> GVCF" | tee -a "$pipeline_log"
        time bash "$script2" \
            "$reference" \
            "$samples" \
            "$bamdir" \
            "$gvcfdir" \
            "$outdir"
    else
        echo "Step 2 already finished, skip." | tee -a "$pipeline_log"
    fi
fi

if [ "$run_step3" = true ]; then
    if [ "$step3_done" = false ]; then
        echo "Step 3: GVCF -> VCF" | tee -a "$pipeline_log"
        time bash "$script3" \
            "$reference" \
            "$samples" \
            "$gvcfdir" \
            "$vcfdir" \
            "$outdir" \
            "$outname"
    else
        echo "Step 3 already finished, skip." | tee -a "$pipeline_log"
    fi
fi

echo "Pipeline command finished." | tee -a "$pipeline_log"
echo "Result summary:" | tee -a "$pipeline_log"
echo "  BAM dir : $bamdir" | tee -a "$pipeline_log"
echo "  GVCF dir: $gvcfdir" | tee -a "$pipeline_log"
echo "  VCF dir : $vcfdir" | tee -a "$pipeline_log"
