wrapper_postprocess <- function(res_link, ancestry_file, params){
  
  # Load res and anc files
  res = unique(read.table(res_link, sep='\t', header=T))
  res = cbind(res[,4:ncol(res)], res[,1:3])
  anc = read.table(ancestry_file, sep='\t', col.names = c('simplesample','ANC'))
  
  print(dim(res))
  print(length(unique(res$start)))
  # Overlap
  resOverlaps = overlap_res_with_cyto_cds_mcnvs_recInvs(res)
  
  # Filer
  resOverlapsFilter = filter_res_t2t_flip(resOverlaps, keep_flip = F)#params$return_full_res)
  
  # If ref is already above cure threshold, we are not interested in further improvements
  resOverlapsFilter[resOverlapsFilter$res_ref > params$cure_threshold_pct,]$mut_max = 'ref'
  #resOverlapsFilter[resOverlapsFilter$res_ref > params$cure_threshold_pct,]$size_mut = 0
  
  # Add a lot of important info columns
  resOverlapsFilterInfo = add_sizemut_invinvolved_etal_to_res(resOverlapsFilter, 
                                                              cure_threshold_pct = params$cure_threshold_pct)
  
  # Add info how many apes are resolved. Once that is done, remove the apes.
  resOverlapsFilterInfoApes = add_n_ape_resolved(resOverlapsFilterInfo, params$ape_samples, res_th = 95)
  resOverlapsFilterInfoApesRemoved = resOverlapsFilterInfoApes[!(resOverlapsFilterInfoApes$sample %in% params$ape_samples),]
  
  # This is the point where names have to shrink.
  rOFIAP = resOverlapsFilterInfoApesRemoved
  

  
  if (params$return_full_res){
    return(rOFIAP)
  }
  # find sSVs
  
  ssvs = res_find_sSVs(rOFIAP, overlap_mode=params$overlap_mode, params$length_limit_bp)
  
  # All members of the ssvs defined earlier
  ssvs_full = rOFIAP[rOFIAP$start %in% unique(ssvs$start),]
  
  # Find those that have eehmm. Invs?
  ssvs_full$overlap =  unlist(lapply(ssvs_full$mut_max, determine_inv_overlapper, overlap_mode = params$overlap_mode))
  
  return(ssvs_full)
  
}

make_inferno_wrapper <- function(ssvs_full_f, params){
  
  print('hi')
  # First plot and the start order
  p1_start_order = plot_overview_n_mut(ssvs_full_f)
  barplot = p1_start_order[[1]]
  start_order = p1_start_order[[2]]
  
  p2_inferno = make_inferno_heatmap(ssvs_full_f, solve_th = params$cure_threshold_pct, process_all = F)
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
  
  print(plot_combine)
  # Consider saving.
  
  #res_locus = plot_one_locus_v2(inferno, anc, start = 1217460, solve_th = 98)
  # res_locus = plot_one_locus_v2(inferno, anc, start = 144039365, solve_th = 98)
  #res_locus = plot_one_locus_v2(inferno, anc, start = 119747586, solve_th = 98)
  if (params$saveInferno){
    ggsave(plot_combine, file = '../plots/Fig2a_raw.pdf', device='pdf', width = 20, height = 10, units='cm')
  }
  
  
}
