#!/bin/bash
set -euo pipefail

echo "argc=$#"
echo "argv=$@"

if [ $# -ne 6 ]; then
    echo "Usage: bash 03_wgs_gvcf_to_vcf.sh <reference.fa> <sample1,sample2,...> <gvcfdir> <vcfdir> <outdir> <outname>"
    echo "Example:"
    echo "bash 03_wgs_gvcf_to_vcf.sh /home/xb/raw_data/Calypte_anna.fna SRR949793,SRR949794 /home/xb/result/gvcf /home/xb/result/vcf /home/xb/result all_samples"
    exit 1
fi

reference="$1"
samples="$2"
gvcfdir="$3"
vcfdir="$4"
outdir="$5"
outname="$6"

logdir="${outdir}/logs/step3"
statusdir="${outdir}/status/step3"
intervaldir="${outdir}/intervals"
partdir="${vcfdir}/parts/${outname}"

mkdir -p "$vcfdir" "$logdir" "$statusdir" "$intervaldir" "$partdir"

if [ ! -f "$reference" ]; then
    echo "ERROR: reference not found: $reference"
    exit 1
fi

if [ ! -d "$gvcfdir" ]; then
    echo "ERROR: gvcfdir not found: $gvcfdir"
    exit 1
fi

if [ ! -f "${reference}.fai" ]; then
    echo "ERROR: fasta index not found: ${reference}.fai"
    exit 1
fi

samples=$(echo "$samples" | tr ',' ' ')

interval_parts="${WGS_INTERVAL_PARTS:-20}"
parallel_jobs="${WGS_PARALLEL_JOBS:-20}"
gatk_java_options="${WGS_GATK_JAVA_OPTIONS:--Xmx24g}"
bcftools_threads="${WGS_BCFTOOLS_THREADS:-8}"
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

if ! [[ "$bcftools_threads" =~ ^[0-9]+$ ]] || [ "$bcftools_threads" -lt 1 ]; then
    echo "ERROR: WGS_BCFTOOLS_THREADS must be a positive integer: $bcftools_threads"
    exit 1
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

log_file="${logdir}/${outname}.log"
done_file="${statusdir}/${outname}.done"

combined_gvcf="${vcfdir}/${outname}.HC.g.vcf.gz"
merged_vcf="${vcfdir}/${outname}.merged.HC.vcf.gz"
snp_vcf="${vcfdir}/${outname}.HC.snp.vcf.gz"
snp_filter_vcf="${vcfdir}/${outname}.HC.snp.filter.vcf.gz"
snp_pass_vcf="${vcfdir}/${outname}.HC.snp.pass.vcf.gz"
invariant_vcf="${vcfdir}/${outname}.invariant.filtered.fixed.clean.vcf.gz"
final_vcf="${vcfdir}/${outname}.novar.noindel.filtered.vcf.gz"
final_tbi="${final_vcf}.tbi"
stats_file="${vcfdir}/${outname}.novar.noindel.filtered.stats.txt"

if [ -f "$final_vcf" ] && [ -f "$final_tbi" ] && [ -f "$done_file" ] && [ -f "$log_file" ] && grep -q "STEP3_SUCCESS" "$log_file"; then
    echo "** final vcf exists and log/status confirmed, skip **"
    exit 0
fi

sample_gvcfs=()
for sample in $samples; do
    gvcf="${gvcfdir}/${sample}.g.vcf.gz"
    gvcf_tbi="${gvcf}.tbi"

    if [ ! -f "$gvcf" ]; then
        echo "ERROR: g.vcf.gz not found: $gvcf"
        exit 1
    fi

    if [ ! -f "$gvcf_tbi" ]; then
        echo "ERROR: g.vcf.gz index not found: $gvcf_tbi"
        exit 1
    fi

    sample_gvcfs+=("-V" "$gvcf")
done

rm -f "$done_file"
rm -f "$combined_gvcf" "${combined_gvcf}.tbi"
rm -f "$merged_vcf" "${merged_vcf}.tbi"
rm -f "$snp_vcf" "${snp_vcf}.tbi"
rm -f "$snp_filter_vcf" "${snp_filter_vcf}.tbi"
rm -f "$snp_pass_vcf" "${snp_pass_vcf}.tbi"
rm -f "$invariant_vcf" "${invariant_vcf}.tbi"
rm -f "$final_vcf" "${final_vcf}.tbi"
rm -f "$stats_file"

{
    echo "[$(date '+%F %T')] outname=$outname start"
    echo "reference=$reference"
    echo "samples=$samples"
    echo "gvcfdir=$gvcfdir"
    echo "vcfdir=$vcfdir"
    echo "interval_mode=$interval_mode"
    echo "interval_parts=$interval_parts"
    echo "parallel_jobs=$parallel_jobs"
    echo "gatk_java_options=$gatk_java_options"
    echo "bcftools_threads=$bcftools_threads"
    echo "interval_files=${interval_files[*]}"

    fail_file="${partdir}/${outname}.failed"
    rm -f "$fail_file"

    part_merged_vcfs=()
    part_no=0
    for interval_file in "${interval_files[@]}"; do
        part_no=$((part_no + 1))
        part_tag=$(printf "part%02d" "$part_no")
        part_gvcf="${partdir}/${outname}.${part_tag}.HC.g.vcf.gz"
        part_vcf="${partdir}/${outname}.${part_tag}.merged.HC.vcf.gz"
        part_log="${logdir}/${outname}.${part_tag}.log"
        part_merged_vcfs+=("$part_vcf")

        if [ -f "$part_vcf" ] && [ -f "${part_vcf}.tbi" ]; then
            echo "[$(date '+%F %T')] ${part_tag} merged VCF exists, skip"
            continue
        fi

        throttle_jobs "$parallel_jobs"
        (
            set -euo pipefail
            echo "[$(date '+%F %T')] ${part_tag} start interval=$interval_file"

            run_gatk CombineGVCFs \
                -R "$reference" \
                -L "$interval_file" \
                "${sample_gvcfs[@]}" \
                -O "$part_gvcf"

            if [ ! -f "${part_gvcf}.tbi" ]; then
                run_gatk IndexFeatureFile -I "$part_gvcf"
            fi

            run_gatk GenotypeGVCFs \
                -all-sites \
                -R "$reference" \
                -V "$part_gvcf" \
                -O "$part_vcf"

            if [ ! -f "${part_vcf}.tbi" ]; then
                run_gatk IndexFeatureFile -I "$part_vcf"
            fi

            echo "[$(date '+%F %T')] ${part_tag} done"
        ) > "$part_log" 2>&1 || {
            echo "${part_tag} failed, see $part_log" >> "$fail_file"
        } &
    done

    if ! wait_for_jobs "$fail_file"; then
        cat "$fail_file"
        exit 1
    fi

    gather_gvcf_args=()
    for interval_file in "${interval_files[@]}"; do
        part_tag=$(basename "$interval_file" .intervals)
        part_tag="${part_tag##*.}"
        part_gvcf="${partdir}/${outname}.${part_tag}.HC.g.vcf.gz"
        if [ -f "$part_gvcf" ]; then
            gather_gvcf_args+=("-I" "$part_gvcf")
        fi
    done

    run_gatk GatherVcfs \
        "${gather_gvcf_args[@]}" \
        -O "$combined_gvcf"

    if [ ! -f "${combined_gvcf}.tbi" ]; then
        run_gatk IndexFeatureFile -I "$combined_gvcf"
    fi

    echo "** ${outname}.HC.g.vcf.gz done **"

    gather_vcf_args=()
    for part_vcf in "${part_merged_vcfs[@]}"; do
        if [ ! -f "$part_vcf" ]; then
            echo "ERROR: missing part VCF: $part_vcf"
            exit 1
        fi
        gather_vcf_args+=("-I" "$part_vcf")
    done

    run_gatk GatherVcfs \
        "${gather_vcf_args[@]}" \
        -O "$merged_vcf"

    if [ ! -f "${merged_vcf}.tbi" ]; then
        run_gatk IndexFeatureFile -I "$merged_vcf"
    fi

    echo "** ${outname}.merged.HC.vcf.gz done **"

    run_gatk SelectVariants \
        -select-type SNP \
        -V "$merged_vcf" \
        -O "$snp_vcf"

    echo "** ${outname}.HC.snp.vcf.gz done **"

    run_gatk VariantFiltration \
        -V "$snp_vcf" \
        --filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
        --filter-name "Filter" \
        -O "$snp_filter_vcf"

    echo "** SNPs filter done **"

    bcftools view --threads "$bcftools_threads" -f PASS "$snp_filter_vcf" -Oz -o "$snp_pass_vcf"
    tabix -p vcf "$snp_pass_vcf"

    echo "** ${outname}.HC.snp.pass.vcf.gz done **"

    bcftools +fill-tags "$merged_vcf" -- -t F_MISSING \
        | bcftools view --threads "$bcftools_threads" -i 'COUNT(GT="alt")=0 && F_MISSING<0.2' \
        | bcftools +fixploidy -- -f 2 \
        | bcftools view --threads "$bcftools_threads" -e 'ALT="*"' \
        -Oz -o "$invariant_vcf"

    tabix -p vcf "$invariant_vcf"

    echo "** ${outname}.invariant.filtered.fixed.clean.vcf.gz done **"

    bcftools concat --threads "$bcftools_threads" -a \
        "$invariant_vcf" \
        "$snp_pass_vcf" \
        -Oz -o "$final_vcf"

    tabix -p vcf "$final_vcf"

    bcftools stats "$final_vcf" > "$stats_file"

    echo "** ${outname}.novar.noindel.filtered.vcf.gz done **"
    echo "STEP3_SUCCESS"
    echo "[$(date '+%F %T')] outname=$outname done"
} > "$log_file" 2>&1

echo "OK $(date '+%F %T')" > "$done_file"
