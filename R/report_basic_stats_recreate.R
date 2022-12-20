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
ape_samples_prim = c("Kamilah_GGO.alt", "Susie_PAB.alt", "Kamilah_GGO.pri", "Susie_PAB.pri","Mhudiblu_PPA.alt",
                     "Mhudiblu_PPA.pri" , "Clint_PTR.alt")
params$ape_samples = c(ape_samples, ape_samples_prim)
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
params$mcnv_mode = 'borders' 

report_enrichment_stats_mcnv(res, unique(ssvs$start), params, middle_only = T)




# Process res to get heatmap.
res2b = add_sizemut_invinvolved_etal_to_res(res, params$cure_threshold_pct)
res2a = res2b
res2a$overlap = unlist(lapply(res2a$mut_max, determine_inv_overlapper, overlap_mode = params$overlap_mode))

res2a[(res2a$overlap == F) & (res2a$result %in% c('2','3')), 'result'] = '1'
res2a[res2a$result == '0', 'result'] = 'Ref'
res2a[res2a$result == '1', 'result'] = '1 SV'
res2a[res2a$result == '2', 'result'] = '2 SVs'
res2a[res2a$result == '3', 'result'] = '3 SVs'
res2a[res2a$result == 'Contig-break', 'result'] = 'No contiguous asm'
res2a[res2a$result == 'Not explained', 'result'] = 'Unexplained'

#res2a$result = as.numeric(res2a$result)

res3 = (res2a[,c('start','sample','result')])
test = melt(res3, is.vars = 'result')
rescast = t(as.matrix(cast(res3, start ~  sample, value.var= result)))
colnames(rescast) = NULL

color_factors = factor(names(table(rescast)), levels = c('No contiguous asm', 'Unexplained','Ref','1 SV', '2 SVs', '3 SVs'))

cols = c('#2d2e30', '#585b61', '#466299', '#8eba7f', '#edb03e', '#ed583e')
colors = structure(cols, names = levels(color_factors))

# Group rows by ancestry

anc = read.table(ancestry_file, sep='\t')
tt = str_split(row.names(rescast), '[.]|_',)
namex = unlist(map(tt, 1))
namex2 = data.frame(V1 = namex)
ancx =  left_join(namex2, anc, by='V1')
row_order = order(ancx$V2)
rescast2 = rescast[row_order,]
ancs =  ancx[row_order,]$V2
# mycolors <- newCols(length(color_factors))
# names(mycolors) <- color_factors
# mycolors <- list(category = mycolors)
#colors = structure(1:length(unique(res3$result)), names = unique(res3$result)) # black, red, green, blue
row.names(rescast2) = paste0(ancs, '___', row.names(rescast2))
Heatmap(rescast2[,sort_order], col = colors, cluster_columns = F, cluster_rows = F,
        column_split = rep(1:5, as.numeric(table(res_slice$category))))

# THIS IS FOR THE 1 ST KIND OF HEATMAP
res2 = add_sizemut_invinvolved_etal_to_res(res, params$cure_threshold_pct)
res2$overlap = unlist(lapply(res2$mut_max, determine_inv_overlapper, overlap_mode = params$overlap_mode))
res2[(res2$overlap == F) & (res2$result %in% c('2','3')), 'result'] = '1'
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
res_slice[((res_slice$n0 + res_slice$n1 + res_slice$n2 + res_slice$n3) == 0) & (res_slice$nbreak >= res_slice$nunexpl), 'category'] = 'No Contig ASM'
res_slice[((res_slice$n0 + res_slice$n1 + res_slice$n2 + res_slice$n3) == 0) & (res_slice$nbreak < res_slice$nunexpl), 'category'] = 'Unexplained SVs'

table(res_slice$category)

res_slice$category = factor(res_slice$category, levels=c('No Contig ASM', 'Unexplained SVs', 'Ref','Simple SVs','sSVs'))
sort_order = order(res_slice$category, res_slice$n0 + res_slice$n1 + res_slice$n2 + res_slice$n3, res_slice$nunexpl, decreasing = F)
res_plot = t(res_slice[sort_order,c('n3', 'n2', 'n1', 'n0',  'nunexpl', 'nbreak')])

white_breaks = cumsum(as.numeric(table(res_slice$category)))
pheatmap(res_plot,
         cluster_cols = F,
         cluster_rows = F,
         gaps_col = white_breaks,
         color = c('#2d2e30', hcl.colors(10, "BluYl")))




