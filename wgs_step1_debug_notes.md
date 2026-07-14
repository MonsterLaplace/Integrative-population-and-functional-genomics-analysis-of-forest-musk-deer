# WGS step1 debug notes

## What changed

- `01_wgs_fastq_to_bam.sh` now writes timestamped `START` and `END exit=...` records around major commands.
- Each sample log now records the exact command that failed, the Bash line number, and the exit code.
- Raw FASTQ files are checked with `gzip -t` before `fastqc` and `fastp`.
- `fastp` thread count is configurable with `FASTP_THREADS`; default is `samtools_threads`.
- `fastp` can be force-stopped with `FASTP_TIMEOUT`; default `0` means no timeout.
- `fastp` minimum length is configurable with `FASTP_MIN_LENGTH`; default is `75`.
- If R1/R2 read length is shorter than `FASTP_MIN_LENGTH`, the script automatically lowers `fastp -l` to the shorter observed read length for that sample and logs a warning.
- `00_run_all_wgs_pipeline.sh` now tees child script output into `logs/pipeline/run_all.log`.
- When one sample fails in step1, the script writes `<outdir>/status/step1/<sample>.failed`, prints the last 40 log lines, and continues with the next sample.
- After all step1 samples finish, the script exits non-zero if any samples failed. This prevents `all` mode from starting step2 with missing BAM files.

## Suggested retry command

Use a timeout while debugging suspected `fastp` hangs. For example, stop `fastp` if it runs longer than 2 hours:

```bash
FASTP_TIMEOUT=7200 FASTP_THREADS=4 bash 00_run_all_wgs_pipeline.sh <datadir> <outdir> <reference.fa> <outname> <aligner_threads> <samtools_threads> 1
```

For short-read datasets, such as 49 bp reads, you can also set the minimum length explicitly:

```bash
FASTP_MIN_LENGTH=30 FASTP_TIMEOUT=7200 FASTP_THREADS=4 bash 00_run_all_wgs_pipeline.sh <datadir> <outdir> <reference.fa> <outname> <aligner_threads> <samtools_threads> 1
```

If a sample hangs or fails, check:

```bash
tail -n 80 <outdir>/logs/step1/<sample>.log
tail -n 120 <outdir>/logs/pipeline/run_all.log
cat <outdir>/status/step1/<sample>.failed
```

Common signals:

- Stops after `START: fastp ...` with no `END`: `fastp` is hanging or the filesystem is stalled.
- `END: timeout ... exit=124`: `fastp` exceeded `FASTP_TIMEOUT`.
- `gzip -t ... exit=1`: the raw FASTQ gzip is corrupt or incomplete.
- `fastp output file missing or empty`: `fastp` exited but did not produce usable clean FASTQ files.
- `too_short_reads` equals nearly all reads in `fastp.json`: `FASTP_MIN_LENGTH` was probably higher than the actual read length.
