###################################################
#1. C、S、M、T score进行窗口评分
###################################################
python integrate_conflict_loci_windows.py \
  --deleterious deleterious_window_density.tsv \
  --selection selection_integrated_windows.tsv \
  --musk-genes musk_function_genes.tsv \
  --roh roh_window_frequency.tsv \
  --tajima tajimaD_50kb.tsv \
  --pixy-pi FMdeer_pixy.merged_pi.txt \
  --annotation FMdeer_unified_gene_annotation.unique.tsv \
  --gene-coords gene_coordinates.tsv \
  -o conflict_loci_integrated_table.tsv

###################################################
#2. 输出窗口关联基因
###################################################

python postprocess_conflict_loci.py \
  --windows conflict_loci_integrated_table.tsv \
  --musk-genes musk_function_genes.tsv \
  --annotation FMdeer_unified_gene_annotation.unique.tsv \
  --gene-coords gene_coordinates.tsv \
  --outdir postprocess_results \
  --top-labels 5

###################################################
#3. 检查各个score窗口值缺失情况，这一步不是必须
###################################################
  for f in C_score.bed S_score.bed M_score.bed T_score.bed
do
  echo "===== $f ====="
  awk '
  NF < 4 {bad++}
  $4 == "NA" || $4 == "NaN" || $4 == "Inf" || $4 == "-Inf" {bad++}
  $4 ~ /^-?[0-9.]+([eE][-+]?[0-9]+)?$/ {
    n++;
    if(n==1 || $4 < min) min=$4;
    if(n==1 || $4 > max) max=$4;
  }
  END {
    print "valid =", n, "bad =", bad+0, "min =", min, "max =", max
  }' $f
done

###################################################
#4. 作图
###################################################
python make_final_quadrant_detail_tables.py
python make_final_quadrant_summary_tables.py
python make_final_top1_tables.py
Rscript plot_conflict_quadrant.R
Rscript plot_quadrant_window_count_pies.R
Rscript plot_window_scores_chromosome_tracks.R
