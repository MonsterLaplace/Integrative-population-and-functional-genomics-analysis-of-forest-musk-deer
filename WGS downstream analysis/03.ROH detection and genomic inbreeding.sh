######################################
#1. ROH计算
######################################
plink --vcf FMdeer.nomaf.recode.vcf.gz \
      --double-id \
      --allow-extra-chr \
      --make-bed \
      --chr-set 33 no-mt \
      --out all_samples.nomaf
wait

plink --bfile all_samples.nomaf \
      --homozyg \
      --homozyg-window-snp 50 \
      --homozyg-snp 50 \
      --homozyg-kb 100 \
      --homozyg-density 50 \
      --homozyg-gap 1000 \
      --homozyg-window-het 1 \
      --allow-extra-chr \
      --chr-set 33 no-mt \
      --out all_samples.roh

##################################
#2. 后续统计处理以及inbreeding计算
##################################
python stat_roh_by_sample.py
python stat_roh_count_and_length.py
python make_roh_window_frequency.py -r  all_samples.roh.hom -f  FMdeer.primary.fa.fai -o roh_window_frequency.tsv
python calc_roh_window_frequency.py
Rscript plot_roh_final.R
