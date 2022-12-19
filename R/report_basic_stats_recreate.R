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

source('report_basic_stats_functions.R')
source('report_basic_stats_wrappers.R')

params = list()

# Define links
res_link = find_reslink(twenty=F)
ancestry_file = '/Users/hoeps/Desktop/desktop_31_july_2022/pics_20_selected/ancestries.tsv'
ape_samples = c('gorGor6', 'panPan3', 'panTro6', 'ponAbe3', 'rheMac10')#,
#                "Clint_PTR.alt", "Mhudiblu_PPA.alt")

params$ape_samples = ape_samples
params$cure_threshold_pct = 98
params$length_limit_bp = 10000
params$make_plots = F
params$overlap_mode='any_overlap'#, 'inv_overlap' #'any_overlap'#
params$return_full_res=T
params$saveInferno = F
res = wrapper_postprocess(res_link, ancestry_file, params)

params$return_full_res=F
ssvs = wrapper_postprocess(res_link, ancestry_file, params)
pp = make_inferno_wrapper(ssvs, params)
# plus25, minus25, plus0, borders
params$mcnv_mode = 'plus25' 
params$mcnv_mode = 'minus25' 
report_enrichment_stats_mcnv(res3, unique(ssvs$start), params, middle_only = T)


# Process res to get heatmap.

res2 = add_sizemut_invinvolved_etal_to_res(res, params$cure_threshold_pct)
res2$overlap = unlist(lapply(res2$mut_max, determine_inv_overlapper, overlap_mode = params$overlap_mode))
res2[res2$overlap == F, 'result'] = '1'
res2[res2$result == '0', 'result'] = 'Ref'
res2[res2$result == '1', 'result'] = '1 SV'
res2[res2$result == '2', 'result'] = '2 SVs'
res2[res2$result == '3', 'result'] = '3 SVs'
res2[res2$result == 'Contig-break', 'result'] = 'No contiguous asm'
res2[res2$result == 'Not explained', 'result'] = 'Unexplained'

res_add = res2 %>% group_by(start) %>% mutate(
  n0 = sum(result=='Ref'),
  n1 = sum(result=='1 SV'),
  n2 = sum(result=='2 SVs'),
  n3 = sum(result=='3 SVs'),
  nbreak = sum(result=='No contiguous asm'),
  nunexpl = sum(result=='Unexplained')
)

res_slice = res_add %>% group_by(start) %>% slice(1)

res_slice$category = 'None'
res_slice[(res_slice$n1 > 0), 'category'] = 'Simple SVs'
res_slice[(res_slice$n2 + res_slice$n3) > 0, 'category'] = 'sSVs'
res_slice[((res_slice$n1 + res_slice$n2 + res_slice$n3) == 0), 'category'] = 'Ref'
res_slice[((res_slice$n0 + res_slice$n1 + res_slice$n2 + res_slice$n3) == 0), 'category'] = 'Unclear'

table(res_slice$category)

res_slice$category = factor(res_slice$category, levels=c('Unclear','Ref','Simple SVs','sSVs'))
sort_order = order(res_slice$category, res_slice$nunexpl, decreasing = F)
res_plot = t(res_slice[sort_order,c('n3', 'n2', 'n1', 'n0',  'nunexpl', 'nbreak')])

white_breaks = cumsum(as.numeric(table(res_slice$category)))
pheatmap(res_plot,
         cluster_cols = F, 
         cluster_rows = F,
         gaps_col = white_breaks)



