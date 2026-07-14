#!/bin/bash
set -euo pipefail

echo "argc=$#"
echo "argv=$@"

if [ $# -ne 5 ]; then
    echo "Usage: bash 02_wgs_bam_to_gvcf.sh <reference.fa> <sample1,sample2,...> <bamdir> <gvcfdir> <outdir>"
    echo "Example:"
    echo "bash 02_wgs_bam_to_gvcf.sh /home/xb/raw_data/Calypte_anna.fna SRR949793,SRR949794 /home/xb/result/bam /home/xb/result/gvcf /home/xb/result"
    exit 1
fi

reference="$1"
samples="$2"
bamdir="$3"
gvcfdir="$4"
outdir="$5"

logdir="${outdir}/logs/step2"
statusdir="${outdir}/status/step2"
intervaldir="${outdir}/intervals"
partdir="${gvcfdir}/parts"

mkdir -p "$gvcfdir" "$logdir" "$statusdir" "$intervaldir" "$partdir"

if [ ! -f "$reference" ]; then
    echo "ERROR: reference not found: $reference"
    exit 1
fi

if [ ! -d "$bamdir" ]; then
    echo "ERROR: bamdir not found: $bamdir"
    exit 1
fi

if [ ! -f "${reference}.fai" ]; then
    echo "ERROR: fasta index not found: ${reference}.fai"
    exit 1
fi

ref_prefix="${reference%.*}"
if [ ! -f "${ref_prefix}.dict" ]; then
    echo "ERROR: sequence dictionary not found: ${ref_prefix}.dict"
    exit 1
fi

samples=$(echo "$samples" | tr ',' ' ')

interval_parts="${WGS_INTERVAL_PARTS:-20}"
parallel_jobs="${WGS_PARALLEL_JOBS:-20}"
hc_threads="${WGS_HC_THREADS:-4}"
gatk_java_options="${WGS_GATK_JAVA_OPTIONS:--Xmx24g}"
interval_mode="${WGS_INTERVAL_MODE:-auto}"

if ! [[ "$interval_parts" =~ ^[0-9]+$ ]] || [ "$interval_parts" -lt 1 ]; then
    echo "ERROR: WGS_INTERVAL_PARTS must be a positive integer: $interval_parts"
    exit 1
fi

case "$interval_mode" in
    auto|balanced|chromosome) ;;
    *)
        echo "ERROR: WGS_INTERVAL_MODE must be auto, balanced, or chromosome: $interval_mode"
        exit 1
        ;;
esac

if [ "$interval_parts" -lt 10 ]; then
    echo "WARN: WGS_INTERVAL_PARTS=$interval_parts is less than recommended 10"
fi

if [ "$interval_parts" -gt 20 ]; then
    echo "WARN: WGS_INTERVAL_PARTS=$interval_parts is greater than recommended 20"
fi

if ! [[ "$parallel_jobs" =~ ^[0-9]+$ ]] || [ "$parallel_jobs" -lt 1 ]; then
    echo "ERROR: WGS_PARALLEL_JOBS must be a positive integer: $parallel_jobs"
    exit 1
fi

if ! [[ "$hc_threads" =~ ^[0-9]+$ ]] || [ "$hc_threads" -lt 1 ]; then
    echo "ERROR: WGS_HC_THREADS must be a positive integer: $hc_threads"
    exit 1
fi

total_hc_threads=$((parallel_jobs * hc_threads))
if [ "$total_hc_threads" -gt 120 ]; then
    echo "WARN: WGS_PARALLEL_JOBS * WGS_HC_THREADS = $total_hc_threads, greater than 120 available threads"
fi

run_gatk() {
    if [ -n "$gatk_java_options" ]; then
        gatk --java-options "$gatk_java_options" "$@"
    else
        gatk "$@"
    fi
}

make_intervals() {
    local fai="$1"
    local n_parts="$2"
    local out_prefix="$3"
    local mode="$4"

    if ls "${out_prefix}".part*.intervals >/dev/null 2>&1; then
        return 0
    fi

    awk -v n="$n_parts" -v prefix="$out_prefix" -v mode="$mode" '
        BEGIN { OFS = "\t" }
        function is_chrom_name(name, lower) {
            lower = tolower(name)
            if (lower ~ /(scaffold|contig|unplaced|unlocalized|random|unknown|unmapped|chr[._-]?un)/) return 0
            if (lower ~ /^chr([0-9]+|x|y|z|w|m|mt)$/) return 1
            if (name ~ /^Chr([0-9]+|[XYxyZWzw]|M|MT)$/) return 1
            if (name ~ /^([0-9]+|[XYxyZWzw]|M|MT)$/) return 1
            return 0
        }
        function interval_line(i, s, e) {
            return chr[i] ":" s "-" e
        }
        {
            chr[NR] = $1
            len[NR] = $2
            total += $2
            chrom[NR] = is_chrom_name($1)
            if (chrom[NR]) {
                chrom_count++
                chrom_total += $2
            }
        }
        END {
            if (NR == 0) exit 1

            resolved_mode = mode
            if (mode == "auto") {
                if (chrom_count > 0 && chrom_count <= 80 && chrom_total / total >= 0.5) {
                    resolved_mode = "chromosome"
                } else {
                    resolved_mode = "balanced"
                }
            }

            summary_file = prefix ".summary.tsv"
            print "mode", resolved_mode > summary_file
            print "total_bp", total >> summary_file
            print "chromosome_candidate_count", chrom_count >> summary_file
            print "chromosome_candidate_bp", chrom_total + 0 >> summary_file

            if (resolved_mode == "chromosome") {
                part = 0
                unplaced_file = ""
                unplaced_bp = 0

                for (i = 1; i <= NR; i++) {
                    if (chrom[i]) {
                        part++
                        file = sprintf("%s.part%02d.intervals", prefix, part)
                        print interval_line(i, 1, len[i]) >> file
                        print sprintf("part%02d", part), chr[i], len[i] >> summary_file
                    }
                }

                for (i = 1; i <= NR; i++) {
                    if (!chrom[i]) {
                        if (unplaced_file == "") {
                            part++
                            unplaced_file = sprintf("%s.part%02d.intervals", prefix, part)
                        }
                        print interval_line(i, 1, len[i]) >> unplaced_file
                        unplaced_bp += len[i]
                    }
                }

                if (unplaced_file != "") {
                    print sprintf("part%02d", part), "unplaced_scaffolds", unplaced_bp >> summary_file
                }
            } else {
                target = int((total + n - 1) / n)
                part = 1
                part_size = 0

                for (i = 1; i <= NR; i++) {
                    pos = 1
                    while (pos <= len[i]) {
                        remaining_in_part = target - part_size
                        if (remaining_in_part < 1) remaining_in_part = target

                        take = len[i] - pos + 1
                        if (take > remaining_in_part && part < n) take = remaining_in_part

                        end = pos + take - 1
                        file = sprintf("%s.part%02d.intervals", prefix, part)
                        print interval_line(i, pos, end) >> file
                        part_bp[part] += take

                        part_size += take
                        pos = end + 1

                        if (part_size >= target && part < n) {
                            part++
                            part_size = 0
                        }
                    }
                }

                for (i = 1; i <= part; i++) {
                    print sprintf("part%02d", i), "balanced", part_bp[i] + 0 >> summary_file
                }
            }
        }
    ' "$fai"
}

throttle_jobs() {
    local max_jobs="$1"
    while [ "$(jobs -rp | wc -l)" -ge "$max_jobs" ]; do
        sleep 5
    done
}

wait_for_jobs() {
    local fail_file="$1"
    wait || true
    if [ -f "$fail_file" ]; then
        return 1
    fi
    return 0
}

interval_prefix="${intervaldir}/genome_${interval_mode}_${interval_parts}"
make_intervals "${reference}.fai" "$interval_parts" "$interval_prefix" "$interval_mode"
mapfile -t interval_files < <(ls "${interval_prefix}".part*.intervals | sort)

if [ "${#interval_files[@]}" -eq 0 ]; then
    echo "ERROR: no interval files created under $intervaldir"
    exit 1
fi

echo "reference=$reference"
echo "samples=$samples"
echo "bamdir=$bamdir"
echo "gvcfdir=$gvcfdir"
echo "outdir=$outdir"
echo "interval_mode=$interval_mode"
echo "interval_parts=$interval_parts"
echo "parallel_jobs=$parallel_jobs"
echo "hc_threads=$hc_threads"
echo "total_hc_threads=$total_hc_threads"
echo "gatk_java_options=$gatk_java_options"
echo "interval_files=${interval_files[*]}"

for sample in $samples; do
    bam="${bamdir}/${sample}.sort.rmdup.bam"
    bai1="${bam}.bai"
    bai2="${bamdir}/${sample}.sort.rmdup.bai"

    gvcf="${gvcfdir}/${sample}.g.vcf.gz"
    gvcf_tbi="${gvcf}.tbi"

    log_file="${logdir}/${sample}.log"
    done_file="${statusdir}/${sample}.done"
    sample_partdir="${partdir}/${sample}"
    mkdir -p "$sample_partdir"

    echo "----------------------------------------"
    echo "sample=$sample"
    echo "bam=$bam"
    echo "gvcf=$gvcf"

    if [ -f "$gvcf" ] && [ -f "$gvcf_tbi" ] && [ -f "$done_file" ] && [ -f "$log_file" ] && grep -q "STEP2_SUCCESS" "$log_file"; then
        echo "** ${sample}.g.vcf.gz exists and log/status confirmed, skip **"
        continue
    fi

    if [ ! -f "$bam" ]; then
        echo "ERROR: bam not found: $bam"
        exit 1
    fi

    if [ ! -f "$bai1" ] && [ ! -f "$bai2" ]; then
        echo "ERROR: bam index not found: $bai1 or $bai2"
        exit 1
    fi

    rm -f "$done_file"
    rm -f "$gvcf" "$gvcf_tbi"

    {
        echo "[$(date '+%F %T')] sample=$sample start"
        echo "reference=$reference"
        echo "bam=$bam"
        echo "gvcf=$gvcf"
        echo "interval_mode=$interval_mode"
        echo "interval_parts=$interval_parts"
        echo "parallel_jobs=$parallel_jobs"
        echo "hc_threads=$hc_threads"
        echo "total_hc_threads=$total_hc_threads"
        echo "gatk_java_options=$gatk_java_options"

        fail_file="${sample_partdir}/${sample}.failed"
        rm -f "$fail_file"

        part_gvcfs=()
        part_no=0
        for interval_file in "${interval_files[@]}"; do
            part_no=$((part_no + 1))
            part_tag=$(printf "part%02d" "$part_no")
            part_gvcf="${sample_partdir}/${sample}.${part_tag}.g.vcf.gz"
            part_log="${logdir}/${sample}.${part_tag}.log"
            part_gvcfs+=("$part_gvcf")

            if [ -f "$part_gvcf" ] && [ -f "${part_gvcf}.tbi" ]; then
                echo "[$(date '+%F %T')] ${sample} ${part_tag} exists, skip"
                continue
            fi

            throttle_jobs "$parallel_jobs"
            (
                set -euo pipefail
                echo "[$(date '+%F %T')] ${sample} ${part_tag} start interval=$interval_file"
                run_gatk HaplotypeCaller \
                    -R "$reference" \
                    -I "$bam" \
                    -L "$interval_file" \
                    -O "$part_gvcf" \
                    -ERC GVCF \
                    --native-pair-hmm-threads "$hc_threads"

                if [ ! -f "${part_gvcf}.tbi" ]; then
                    run_gatk IndexFeatureFile -I "$part_gvcf"
                fi

                echo "[$(date '+%F %T')] ${sample} ${part_tag} done"
            ) > "$part_log" 2>&1 || {
                echo "${sample} ${part_tag} failed, see $part_log" >> "$fail_file"
            } &
        done

        if ! wait_for_jobs "$fail_file"; then
            cat "$fail_file"
            exit 1
        fi

        gather_args=()
        for part_gvcf in "${part_gvcfs[@]}"; do
            if [ ! -f "$part_gvcf" ]; then
                echo "ERROR: missing part gVCF: $part_gvcf"
                exit 1
            fi
            gather_args+=("-I" "$part_gvcf")
        done

        run_gatk GatherVcfs \
            "${gather_args[@]}" \
            -O "$gvcf"

        if [ ! -f "$gvcf_tbi" ]; then
            run_gatk IndexFeatureFile -I "$gvcf"
        fi

        echo "STEP2_SUCCESS"
        echo "[$(date '+%F %T')] sample=$sample done"
    } > "$log_file" 2>&1

    echo "OK $(date '+%F %T')" > "$done_file"

    echo "** ${sample}.g.vcf.gz done **"
done
