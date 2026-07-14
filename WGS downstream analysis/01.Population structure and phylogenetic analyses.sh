#1. 群体层面进一步过滤
vcftools \
    --gzvcf $project/02.WGS/06.filter/all_samples.filtered.snp.PASS.vcf.gz \
    --max-missing 0.9 \
    --maf 0.05 \
    --minQ 30 \
    --recode --recode-INFO-all \
    --out $project/02.WGS/06.filter/all_samples.strict
bgzip -f $project/02.WGS/06.filter/all_samples.strict.recode.vcf
tabix -p vcf $project/02.WGS/06.filter/all_samples.strict.recode.vcf.gz

#2.转PLINK格式
plink --vcf $project/02.WGS/06.filter/all_samples.strict.recode.vcf.gz \
      --double-id \
      --allow-extra-chr \
      --make-bed \
      --out $project/02.WGS/07.population/all_samples.strict

#3. LD pruning
plink --bfile $project/02.WGS/07.population/all_samples.strict \
      --indep-pairwise 50 10 0.2 \
      --allow-extra-chr \
      --out $project/02.WGS/07.population/all_samples.prune

plink --bfile $project/02.WGS/07.population/all_samples.strict \
      --extract $project/02.WGS/07.population/all_samples.prune.prune.in \
      --make-bed \
      --allow-extra-chr \
      --out $project/02.WGS/07.population/all_samples.pruned

#4. PCA
VCF2PCACluster -InVCF /data/xb/FMdeer/04.finalSNP/FMdeer.strict.recode.vcf.gz -OutPut FMdeer -InSampleGroup pop.info -Threads 80

#5. ADMIXTURE
for K in 2 3 4 5 6; do
    admixture --cv $project/02.WGS/07.population/all_samples.pruned.bed $K | tee $project/02.WGS/07.population/log${K}.out
done

#5.  构建系统树
python vcf2phylip.py -i $project/02.WGS/06.filter/all_samples.strict.recode.vcf.gz

iqtree2 -s all_samples.strict.recode.min4.phy -m MFP -bb 1000 -nt AUTO -st DNA


