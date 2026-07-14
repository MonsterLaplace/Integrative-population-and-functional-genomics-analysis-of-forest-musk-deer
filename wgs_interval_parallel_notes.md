# WGS interval parallel notes

## What changed

- Step2 (`02_wgs_bam_to_gvcf.sh`) now splits the reference genome into interval files using `<reference.fa>.fai`.
- Interval splitting uses `WGS_INTERVAL_MODE=auto` by default:
  - chromosome-level assembly: chromosome-like contigs are split one per interval file; unplaced/unlocalized/random/scaffold contigs are merged into one final interval file.
  - scaffold-level assembly: contigs are split by cumulative bp size into similarly sized interval files, preserving reference order for safe VCF gathering.
- Step2 runs `gatk HaplotypeCaller` by sample and interval in parallel, then gathers parts back to the original output:
  - `<outdir>/gvcf/<sample>.g.vcf.gz`
- Step3 (`03_wgs_gvcf_to_vcf.sh`) runs `gatk CombineGVCFs` and `gatk GenotypeGVCFs` by interval in parallel, then gathers parts back to:
  - `<outdir>/vcf/<outname>.HC.g.vcf.gz`
  - `<outdir>/vcf/<outname>.merged.HC.vcf.gz`
- Downstream SNP filtering, invariant-site filtering, final concat, tabix, and stats keep the original output names.

## Runtime knobs

Recommended range for genome chunks: 10-20.

For your 120-thread / 1 TB RAM node, the new defaults are:

- `WGS_INTERVAL_PARTS=20`
- `WGS_INTERVAL_MODE=auto`
- `WGS_PARALLEL_JOBS=20`
- `WGS_HC_THREADS=4`
- `WGS_GATK_JAVA_OPTIONS="-Xmx24g"`
- `WGS_BCFTOOLS_THREADS=8`

```bash
WGS_INTERVAL_MODE=auto WGS_INTERVAL_PARTS=20 WGS_PARALLEL_JOBS=20 WGS_HC_THREADS=4 WGS_GATK_JAVA_OPTIONS="-Xmx24g" WGS_BCFTOOLS_THREADS=8 bash 00_run_all_wgs_pipeline.sh \
  <datadir> <outdir> <reference.fa> <outname> <aligner_threads> <samtools_threads> 2,3
```

- `WGS_INTERVAL_MODE`: `auto`, `balanced`, or `chromosome`. Default: `auto`.
- `WGS_INTERVAL_PARTS`: target number of balanced genome chunks. Default: `20`.
- `WGS_PARALLEL_JOBS`: maximum interval jobs running at once. Default: `20`.
- `WGS_HC_THREADS`: native PairHMM threads per HaplotypeCaller job. Default: `4`.
- `WGS_GATK_JAVA_OPTIONS`: Java heap/options passed to GATK. Default: `-Xmx24g`.
- `WGS_BCFTOOLS_THREADS`: compression/filtering threads for bcftools in step3. Default: `8`.

This uses about 80 HaplotypeCaller compute threads in step2 (`20 x 4`) and leaves CPU headroom for Java, disk IO, compression, and system overhead.

If memory is tight, reduce `WGS_PARALLEL_JOBS` first. If CPU is idle and memory is enough, increase `WGS_HC_THREADS` to `5` or `6`, but keep `WGS_PARALLEL_JOBS x WGS_HC_THREADS` near or below 120.

Conservative run:

```bash
WGS_INTERVAL_MODE=auto WGS_INTERVAL_PARTS=20 WGS_PARALLEL_JOBS=12 WGS_HC_THREADS=4 WGS_GATK_JAVA_OPTIONS="-Xmx24g" bash 00_run_all_wgs_pipeline.sh \
  <datadir> <outdir> <reference.fa> <outname> <aligner_threads> <samtools_threads> 2,3
```

More aggressive run:

```bash
WGS_INTERVAL_MODE=auto WGS_INTERVAL_PARTS=20 WGS_PARALLEL_JOBS=20 WGS_HC_THREADS=6 WGS_GATK_JAVA_OPTIONS="-Xmx24g" bash 00_run_all_wgs_pipeline.sh \
  <datadir> <outdir> <reference.fa> <outname> <aligner_threads> <samtools_threads> 2,3
```

## New intermediate locations

- Intervals:
  - `<outdir>/intervals/genome_<mode>_<N>.partXX.intervals`
  - `<outdir>/intervals/genome_<mode>_<N>.summary.tsv`
- Step2 part gVCFs:
  - `<outdir>/gvcf/parts/<sample>/<sample>.partXX.g.vcf.gz`
- Step3 part VCFs:
  - `<outdir>/vcf/parts/<outname>/<outname>.partXX.*.vcf.gz`

Existing completion checks still use the original final files and `STEP2_SUCCESS` / `STEP3_SUCCESS`.
