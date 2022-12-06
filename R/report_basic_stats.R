# Whoeps, 6th Apr 2022
# Explore output file.

library(dplyr)
library(stringr)
library(ggplot2)
library(tidyr)
library(cowplot)

source('report_basic_stats_functions.R')



res_link = find_reslink(twenty=F)
ancestry_file = '/Users/hoeps/Desktop/desktop_31_july_2022/pics_20_selected/ancestries.tsv'

ape_samples = c('gorGor6_hg38', 'panPan3_hg38', 'panTro6_hg38', 'ponAbe3_hg38', 'rheMac10_hg38')
cure_threshold_pct = 98
length_limit_bp = 10000
make_plots = F
overlap_mode='inv_overlap'

# Start computing.
res = unique(read.table(res_link, sep='\t', header=T))

res = overlap_with_cyto_bands(res)
res = overlap_with_core_dups(res)
res = overlap_with_recurrent_invs(res)
res = overlap_with_mcnvs(res)

anc = read.table(ancestry_file, sep='\t', col.names = c('simplesample','ANC'))

res = res[res$end - res$start > length_limit_bp,]

res = res[res$flip_unsure == F,]

res$mut_max = unlist(lapply(res$mut_max, filter_t))

res_plus = enrich_res_with_info(res, cure_threshold_pct = cure_threshold_pct)

# Remove that weird sample
res_plus = res_plus[res_plus$sample != 'HG02666_chrY_hg38',]

#
res_plus = add_n_ape_resolved(res_plus, res_th = 95)

# Filter #2: remove apes
res_plus = res_plus[!(res_plus$sample %in% ape_samples),]

# Drastic Eingriff here
res_plus[res_plus$res_ref > cure_threshold_pct,]$mut_max = 'ref'
res_plus[res_plus$res_ref > cure_threshold_pct,]$size_mut = 0

r2 = res_plus[res_plus$start==144039365,]
head(r2)

if (make_plots){
  # Add info


  p_all = plot_overview_n_mut(res_plus, sort_by_ref = F)

  plot_one_locus_v1(res_plus, start = 4186831)
}




nsvs = res_find_nSVs(res_plus, overlap_mode=overlap_mode)


# Experiment

evaluate(res, nsvs)

nsvs_full = res_plus[res_plus$start %in% unique(nsvs$start),]
nsvs_full$overlap =  unlist(lapply(nsvs_full$mut_max, determine_inv_overlapper, overlap_mode = overlap_mode))

p1_start_order = plot_overview_n_mut(nsvs_full)
p1 = p1_start_order[[1]]
start_order = p1_start_order[[2]]
plot(1,1)

p2_inferno = make_inferno_heatmap(nsvs_full, solve_th = cure_threshold_pct, process_all = T)
p2 = p2_inferno[[1]]
inferno = p2_inferno[[2]]




res_locus = plot_one_locus_v2(inferno, anc, start = 1217460, solve_th = 98)
res_locus = plot_one_locus_v2(inferno, anc, start = 144039365, solve_th = 98)
res_locus = plot_one_locus_v2(inferno, anc, start = 119747586, solve_th = 98)
res_locus = plot_one_locus_v2(inferno, anc, start = 4186831, solve_th = 98)

# for (start in unique(inferno$start)){
#
#   outplot = plot_one_locus_v2(inferno, anc, start = start, solve_th = 98)[[1]]
#   save_plot(outplot, filename = paste('../diagramplots/',start, '.pdf'), device='pdf', base_height = 5, base_width = 7)
#
#   print('sup')
# }














if (F){
  library(cowplot)
  save_plot(p1, filename = '2a_bars_2.pdf', device='pdf', base_height = 10, base_width = 7)
  save_plot(p2, filename = '2b.pdf_2.pdf', device='pdf', base_height = 10, base_width = 12)

}
if (F){
plot_overview_n_mut(res_plus[res_plus$start %in% unique(nsvs$start),])
#uniqs = nsvs %>% group_by(start) %>% slice(1)
#uniqs = uniqs[(uniqs$end - uniqs$start) > length_limit_bp,]

# Filter #8: branch out from #7

# Filter #9: branch from #8



# Filter #10: branch from $9
a = uniqs[c('seqname','start','end')]

# Filter #11
res_nobreak = res[(res$exceeds_x == F) &  (res$exceeds_y == F),]


# Make x-order
# res = (res %>% group_by(start) %>% mutate(missing = n_asms - length(result)))



}





