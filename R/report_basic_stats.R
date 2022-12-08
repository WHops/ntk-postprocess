# Whoeps, 6th Apr 2022
# Explore output file.

library(dplyr)
library(stringr)
library(ggplot2)
library(tidyr)
library(cowplot)
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(patchwork)
library(ggplotify)
library(reshape2)
library(dplyr)
source('report_basic_stats_functions.R')

ssv_starts = c( 135571499, 154106432 , 49223962,  56758523 , 17939968,
                20998524 , 54689990 ,  7009751 , 35983478 , 37739530 ,
                43203350 , 14200886  ,21003474  ,28206731 , 72044607 ,
                82231548 , 84256360  ,18581382 ,  4186831 , 89481973 ,
                45338205 , 87850116 , 37623508  , 5735431 ,102392422 ,
                143470184  ,72293822 , 17456276 ,195515228 ,103543499 ,
                108148088 ,109662781, 119747586,12830558 ,207485628,
                20991472)

# Define links
res_link = find_reslink(twenty=F)
ancestry_file = '/Users/hoeps/Desktop/desktop_31_july_2022/pics_20_selected/ancestries.tsv'

# Hardcoded parameters
#ape_samples = c('gorGor6_hg38', 'panPan3_hg38', 'panTro6_hg38', 'ponAbe3_hg38', 'rheMac10_hg38')
ape_samples = c('gorGor6', 'panPan3', 'panTro6', 'ponAbe3', 'rheMac10')

cure_threshold_pct = 98
length_limit_bp = 10000
make_plots = F
overlap_mode='any_overlap'#, 'inv_overlap' #'any_overlap'#


# Start computing.
res = unique(read.table(res_link, sep='\t', header=T))
res = cbind(res[,4:ncol(res)], res[,1:3])
anc = read.table(ancestry_file, sep='\t', col.names = c('simplesample','ANC'))

# One more try. What if we exclude everything that has not got any resolution? 
res = res[!(res$sample %in% ape_samples),]
res = res[res$sample != 'T2T',]
res2 =res %>% group_by(start) %>% mutate(n_fixed = (sum(res_max > cure_threshold_pct)))
res3 = res[res2$n_fixed > 0,]

# Stuff that got fixed by at least one non-ref SV. 
res_nonref = res %>% group_by(start) %>% mutate(n_nonref = (sum((res_max > cure_threshold_pct) & (mut_max != 'ref')) ))
res4 = res[res_nonref$n_nonref == 0,]

# Stuff that has no contig break in at least 50% of cases.
res_nonref = res %>% group_by(start) %>% mutate(n_exceeds = (sum( exceeds_y)))
res_nonref2 = res[res_nonref$n_exceeds < 55,]

mcnv_link = '/Users/hoeps/PhD/projects/nahrcall/analyses_paper/mcnv_enrichment/data/all_mcnvs_merged_plus25pct.bed'
mcnv_link = '/Users/hoeps/PhD/projects/nahrcall/analyses_paper/mcnv_enrichment/data/all_mcnvs_merged_minus25pct.bed'
mcnv_link = '/Users/hoeps/PhD/projects/nahrcall/analyses_paper/mcnv_enrichment/data/all_mcnvs_merged.bed'

mcnv_link = '/Users/hoeps/PhD/projects/nahrcall/analyses_paper/mcnv_enrichment/data/borders_only.bed'
report_enrichment_stats_mcnv(res3, ssv_starts, mcnv_link, middle_only = T)

####### All for pre-processing the data in some ways ###########

# Overlap res with other data.
res = overlap_with_cyto_bands(res)
res = overlap_with_core_dups(res)
res = overlap_with_recurrent_invs(res)
res = overlap_with_mcnvs(res)




# Filter res
res = res[res$end - res$start > length_limit_bp,]
res = res[res$flip_unsure == F,]

# Ummmm
res$mut_max = unlist(lapply(res$mut_max, filter_t))
res_plus = enrich_res_with_info(res, cure_threshold_pct = cure_threshold_pct)

# Remove that weird sample
res_plus = res_plus[res_plus$sample != 'HG02666_chrY_hg38',]

# Add info how many apes are resolved. Once that is done, remove the apes.
res_plus = add_n_ape_resolved(res_plus, res_th = 95)
#res_plus$n_ape_resolved = 5
res_plus = res_plus[!(res_plus$sample %in% ape_samples),]

# If ref is already above cure threshold, we are not interested in further improvements
res_plus[res_plus$res_ref > cure_threshold_pct,]$mut_max = 'ref'
res_plus[res_plus$res_ref > cure_threshold_pct,]$size_mut = 0

# find sSVs
ssvs = res_find_sSVs(res_plus, overlap_mode=overlap_mode)

####### All for pre-processing the data in some ways ###########

# This function checks if the ssvs that we expect from visual inspection are here.
evaluate(res, ssvs)

# All members of the ssvs defined earlier
ssvs_full = res_plus[res_plus$start %in% unique(ssvs$start),]

# Find those that have eehmm. Invs?
ssvs_full$overlap =  unlist(lapply(ssvs_full$mut_max, determine_inv_overlapper, overlap_mode = overlap_mode))

# First plot and the start order
p1_start_order = plot_overview_n_mut(ssvs_full)
barplot = p1_start_order[[1]]
start_order = p1_start_order[[2]]

p2_inferno = make_inferno_heatmap(ssvs_full, solve_th = cure_threshold_pct, process_all = F)
heatmap = p2_inferno[[1]]
inferno = p2_inferno[[2]]

plot.1 = barplot
plot.2 = as.ggplot(heatmap)

layout = "
AABBBBBBBBBBBB
AABBBBBBBBBBBB
AABBBBBBBBBBBB
AABBBBBBBBBBBB
AABBBBBBBBBBBB
AABBBBBBBBBBBB
AABBBBBBBBBBBB
AABBBBBBBBBBBB
AABBBBBBBBBBBB
AABBBBBBBBBBBB
AABBBBBBBBBBBB
AABBBBBBBBBBBB
AABBBBBBBBBBBB
AABBBBBBBBBBBB
AABBBBBBBBBBBB
AABBBBBBBBBBBB
##BBBBBBBBBBBB
##BBBBBBBBBBBB
"


plot_combine = plot.1 + theme(panel.border = element_blank()) +
  plot.2 + plot_layout(guides = "collect", design=layout)#widths = c(0.1,0.9), heights=c(2,1))
plot_combine
# Consider saving.
#ggsave(plot_combine, file = '../plots/Fig2a_raw.pdf', device='pdf', width = 20, height = 10, units='cm')

#res_locus = plot_one_locus_v2(inferno, anc, start = 1217460, solve_th = 98)
# res_locus = plot_one_locus_v2(inferno, anc, start = 144039365, solve_th = 98)
#res_locus = plot_one_locus_v2(inferno, anc, start = 119747586, solve_th = 98)
#res_locus = plot_one_locus_v2(inferno, anc, start = 4186831, solve_th = 98)
save_plot(outplot, filename = paste('../diagramplots/',start, '.pdf'), device='pdf', base_height = 5, base_width = 7)
# for (start in unique(inferno$start)){
#
#   outplot = plot_one_locus_v2(inferno, anc, start = start, solve_th = 98)[[1]]
#   save_plot(outplot, filename = paste('../diagramplots/',start, '.pdf'), device='pdf', base_height = 5, base_width = 7)
#
#   print('sup')
# }

#ssv_folders_to_get = unique(paste(ssvs_full[,c('')]))

link_to_remote_res = '~/cluster21/g/korbel/hoeps/projects/nahr/ntk_interesting_loci/ntk-scan-snakemake/res'
link_to_local_target = '/Users/hoeps/PhD/projects/nahrcall/analyses_paper_2/first_run/res'
out_script = '/Users/hoeps/PhD/projects/nahrcall/analyses_paper_2/first_run/copy_command_auto.sh'

copy_ssvs_pdf_to_local(ssvs_full, link_to_remote_res, link_to_local_target, out_script, execute = F)


# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# if (F){
#   library(cowplot)
#   save_plot(p1, filename = '2a_bars_2.pdf', device='pdf', base_height = 10, base_width = 7)
#   save_plot(p2, filename = '2b.pdf_2.pdf', device='pdf', base_height = 10, base_width = 12)
# 
# }
# if (F){
# plot_overview_n_mut(res_plus[res_plus$start %in% unique(ssvs$start),])
# #uniqs = ssvs %>% group_by(start) %>% slice(1)
# #uniqs = uniqs[(uniqs$end - uniqs$start) > length_limit_bp,]
# 
# # Filter #8: branch out from #7
# 
# # Filter #9: branch from #8
# 
# 
# 
# # Filter #10: branch from $9
# a = uniqs[c('seqname','start','end')]
# 
# # Filter #11
# res_nobreak = res[(res$exceeds_x == F) &  (res$exceeds_y == F),]
# 
# 
# # Make x-order
# # res = (res %>% group_by(start) %>% mutate(missing = n_asms - length(result)))
# 
# 
# 
# }



# Unused

# if (make_plots){
#   # Barplot summarizing everything
#   p_bars = plot_overview_n_mut(res_plus, sort_by_ref = F)
#
#   # Individual plot
#   plot_one_locus_v1(res_plus, start = 4186831)
# }


