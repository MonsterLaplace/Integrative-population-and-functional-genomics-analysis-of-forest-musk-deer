########################
#1. genetic load
########################
java -jar /data/xb/FMdeer/04.finalSNP/geneticload/snpEff/snpEff.jar build -gff3 -v FMdeer -d -noCheckCds -noCheckProtein

java -Xmx20g -jar snpEff.jar -v FMdeer \
    $project/02.WGS/06.filter/all_samples.filtered.snp.PASS.vcf.gz \
    > $project/02.WGS/09.load/all_samples.snpeff.ann.vcf

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/ANN\n' \
  $project/02.WGS/09.load/all_samples.snpeff.ann.vcf.gz \
  > $project/02.WGS/09.load/all_samples.ANN.tsv

########################
#2. 下游分析
########################
#（1）提取各个遗传负荷子集
python extract_functional_variants_with_grantham.py \
  -i all_samples.ANN.tsv \
  -o functional_variants_with_grantham.tsv \
  --severe-threshold 150

awk -F'\t' 'NR>1 && $19=="LoF"{print $1"\t"$2}' functional_variants_with_grantham.tsv > LoF.sites

awk -F'\t' 'NR>1 && $19=="missense_all"{print $1"\t"$2}' functional_variants_with_grantham.tsv > missense_all.sites

awk -F'\t' 'NR>1 && $19=="missense_severe"{print $1"\t"$2}' functional_variants_with_grantham.tsv > missense_severe.sites

awk -F'\t' 'NR>1 && $19=="synonymous"{print $1"\t"$2}' functional_variants_with_grantham.tsv > synonymous.sites

bcftools view -R LoF.sites all_samples.snpeff.ann.vcf.gz -Oz -o LoF.vcf.gz
tabix -p vcf LoF.vcf.gz

bcftools view -R missense_all.sites all_samples.snpeff.ann.vcf.gz -Oz -o missense_all.vcf.gz
tabix -p vcf missense_all.vcf.gz

bcftools view -R missense_severe.sites all_samples.snpeff.ann.vcf.gz -Oz -o missense_severe.vcf.gz
tabix -p vcf missense_severe.vcf.gz

bcftools view -R synonymous.sites all_samples.snpeff.ann.vcf.gz -Oz -o synonymous.vcf.gz
tabix -p vcf synonymous.vcf.gz

#（2）位点级去重
python collapse_functional_variants_to_site_level.py \
  -i functional_variants_with_grantham.tsv \
  -o functional_variants_site_level.tsv

python - <<'PY'
import pandas as pd
df = pd.read_csv("functional_variants_site_level.tsv", sep="\t", dtype=str)

for cls in ["LoF", "missense_severe", "missense_all", "synonymous"]:
    sub = df[df["class"] == cls][["CHROM","POS"]].drop_duplicates()
    sub.to_csv(f"{cls}.site_level.sites", sep="\t", index=False, header=False)
PY

# （3）计算遗传负荷
python calculate_burden_from_site_level_and_vcf.py \
  -s functional_variants_site_level.tsv \
  -v all_samples.snpeff.ann.vcf.gz \
  -o genetic_load_per_sample.tsv \
  --classes LoF missense_severe missense_all

# （4）下游分析与作图（和ROH关联分析）
Rscript compare_genetic_load_domestic_vs_wild.R

python genetic_load_inside_vs_outside_roh.py \
  -s functional_variants_site_level.tsv \
  -v all_samples.snpeff.ann.vcf.gz \
  -r all_samples.roh.bed \
  -o genetic_load_inside_outside_ROH.tsv \
  --classes LoF missense_severe missense_all

Rscript plot_genetic_load_inside_outside_roh.R
Rscript plot_deleterious_density_genomewide.R


  
  
