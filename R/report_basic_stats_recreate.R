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

res = wrapper_postprocess(res_link, ancestry_file, params)

params$return_full_res=F
ssvs = wrapper_postprocess(res_link, ancestry_file, params)

# plus25, minus25, plus0, borders
params$mcnv_mode = 'plus25' 
params$mcnv_mode = 'minus25' 
report_enrichment_stats_mcnv(res3, unique(ssvs$start), params, middle_only = T)
