#####################################
whole WGS pipeline written by xubo
#####################################

使用前请仔细阅读：
该脚本我做了工程化封装，目前分为三大块：1：从fastq文件到bam文件，2：从bam到g.vcf，3：从g.vcf到两套vcf (包含变异和非变异位点)
脚本支持断点续跑，支持跑其中任意一步：

推荐这样组织：
1）脚本目录
/home/xb/scripts/
├── 00_run_all_wgs_pipeline.sh
├── 01_wgs_fastq_to_bam.sh
├── 02_wgs_bam_to_gvcf.sh
└── 03_wgs_gvcf_to_vcf.sh

2）原始数据目录
/home/xb/raw_data/
├── Calypte_anna.fna
├── SRR949793_1.fastq.gz
├── SRR949793_2.fastq.gz
├── SRR949794_1.fastq.gz
└── SRR949794_2.fastq.gz

3）结果目录
/home/xb/result/
├── bam/
├── gvcf/
├── vcf/
├── logs/
│   ├── pipeline/
│   ├── step1/
│   ├── step2/
│   └── step3/
└── status/
    ├── step1/
    ├── step2/
    └── step3/

使用示例：

可以用：bash /home/xb/scripts/00_run_all_wgs_pipeline.sh -h，或者bash /home/xb/scripts/00_run_all_wgs_pipeline.sh --help查看帮助

跑完整流程
bash /home/xb/scripts/00_run_all_wgs_pipeline.sh \
/home/xb/raw_data \
/home/xb/result \
Calypte_anna.fna \
all_samples \
16 \
4 \
all

这个脚本意思是：
调用总控脚本 00_run_all_wgs_pipeline.sh，
以 /home/xb/raw_data 里的参考基因组和 FASTQ 数据为输入，
把结果输出到 /home/xb/result，
用 Calypte_anna.fna 作为参考基因组，
最终输出文件前缀命名为 all_samples，
bwa-mem2 用 16 线程，samtools 用 4 线程，
并运行完整流程 all（即 step1 + step2 + step3）。

如果只跑第 1 步
bash /home/xb/scripts/00_run_all_wgs_pipeline.sh \
/home/xb/raw_data \
/home/xb/result \
Calypte_anna.fna \
all_samples \
16 \
4 \
1

从第 2 步开始：
bash /home/xb/scripts/00_run_all_wgs_pipeline.sh \
/home/xb/raw_data \
/home/xb/result \
Calypte_anna.fna \
all_samples \
16 \
4 \
2,3

关于断点续跑的问题：
1.统一完成判断规则吗，通过检测文件是否存在，以及日志这一步是否完成来判断文件完整性
用“完成标记文件 + 独立日志文件”双保险
对于每一步/每个样本：
开始运行时写一个日志
成功结束后写一个 done 标记文件
下次判断是否跳过时，不只看结果文件，还要看：
结果文件存在
索引文件存在
done 标记文件存在
done 标记文件中有明确成功标志

日志用来排错
done 标记用来判断完成
例如：

BASH
outdir/logs/step1/sample.step1.log
outdir/status/step1/sample.step1.done

done 文件内容如下：
BASH
OK
2026-06-05 12:30:11
这样判断时非常简单可靠。

另外：
Step1：FASTQ -> BAM
某样本完成条件：
bam/${sample}.sort.rmdup.bam 存在
bam/${sample}.sort.rmdup.bam.bai 存在
status/step1/${sample}.done 存在
logs/step1/${sample}.log 中包含：
BASH
STEP1_SUCCESS

Step2：BAM -> GVCF
某样本完成条件：
gvcf/${sample}.g.vcf.gz 存在
gvcf/${sample}.g.vcf.gz.tbi 存在
status/step2/${sample}.done 存在
logs/step2/${sample}.log 中包含：
BASH
STEP2_SUCCESS

Step3：GVCF -> VCF
整体完成条件：
vcf/${outname}.novar.noindel.filtered.vcf.gz 存在
vcf/${outname}.novar.noindel.filtered.vcf.gz.tbi 存在
status/step3/${outname}.done 存在
logs/step3/${outname}.log 中包含：
BASH
STEP3_SUCCESS

