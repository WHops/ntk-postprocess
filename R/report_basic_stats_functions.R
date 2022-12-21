filter_res_t2t_flip <- function(res_f, ape_samples, keep_flip){
  # Filter res: kick out apes, T2T, and flip_unsures
  #res_f = res_f[!(res_f$sample %in% ape_samples),]
  res_f = res_f[res_f$sample != 'T2T-CHM13v2.0',]
  
  if (!(keep_flip)){
    res_f = res_f[res_f$flip_unsure == F,]
  }
  # Remove that weird sample in case it's there...
  res_f = res_f[res_f$sample != 'HG02666_chrY_hg38',]
  
  return(res_f)
}



overlap_res_with_cyto_cds_mcnvs_recInvs <- function(res_f){
  # Overlap res with other data.
  res_f = overlap_with_cyto_bands(res_f)
  res_f = overlap_with_core_dups(res_f)
  res_f = overlap_with_recurrent_invs(res_f)
  res_f = overlap_with_mcnvs(res_f)
  
  return(res_f)
}

report_enrichment_stats_mcnv <- function(res_f, ssv_starts, params, middle_only = F){
  
  res = res_f[!(res_f$sample %in% params$ape_samples),]
  res = res[res$sample != 'T2T',]
  res2 =res %>% group_by(start) %>% mutate(n_fixed = (sum(res_max > params$cure_threshold_pct)))
  res3 = res[res2$n_fixed > 0,]
  
  if (params$mcnv_mode == 'plus25'){
    mcnv_link = '/Users/hoeps/PhD/projects/nahrcall/analyses_paper/mcnv_enrichment/data/all_mcnvs_merged_plus25pct.bed'
  } else if (params$mcnv_mode == 'minus25'){
    mcnv_link = '/Users/hoeps/PhD/projects/nahrcall/analyses_paper/mcnv_enrichment/data/all_mcnvs_merged_minus25pct.bed'
  } else if (params$mcnv_mode == 'plus0'){
    mcnv_link = '/Users/hoeps/PhD/projects/nahrcall/analyses_paper/mcnv_enrichment/data/all_mcnvs_merged.bed'
  } else if (params$mcnv_mode == 'borders'){
    mcnv_link = '/Users/hoeps/PhD/projects/nahrcall/analyses_paper/mcnv_enrichment/data/borders_only.bed'
  } else {
    break()
  }
  
  
  library(dplyr)
  overlap_tmp_link = 'mcnv_overlappers.bed'
  bedtools_link = '/usr/local/bin/bedtools'
  res_tmp_link = 'res_tmp.tmp'
  
  res_f$oldstart = res_f$start
  if (middle_only == T){
    res_f$start = as.integer(res_f$start + ((res_f$end - res_f$start) / 2))
    res_f$end = res_f$start + 1
  }
  
  # Save res table, so bedtools can work with it.
  write.table(res_f[,1:3], res_tmp_link, sep='\t', col.names = F, row.names = F, quote = F)
  
  # Run bedtools to find intersections betwres_locus_meltn res and coredups
  bed_command = paste0(bedtools_link ,' intersect -wa -a ', res_tmp_link, ' -b ', mcnv_link, ' | sort | uniq > ', overlap_tmp_link)
  system(bed_command)
  
  # Load the results.
  mcnv_overlapping_segments = read.table(overlap_tmp_link, sep='\t')
  colnames(mcnv_overlapping_segments) = c('seqname','start','end')
  mcnv_overlapping_segments$mcnv = T
  
  # Add this information to the original res directory
  res_f$mcnv = NULL
  res_f = dplyr::left_join(res_f, mcnv_overlapping_segments, by=c('seqname', 'start','end'))
  res_f[is.na(res_f$mcnv),'mcnv'] = F
  
  # First, normal:
  res_uniq = res_f %>% group_by(start) %>% slice(1)
  table1 = table(res_uniq$mcnv)
  print(table(res_uniq$mcnv))
  
  # Second, ssvs:
  res_ssvs = res_f[as.numeric(res_f$oldstart) %in% ssv_starts,]
  res_ssvs_uniq = res_ssvs %>% group_by(start) %>% slice(1)
  table2 = table(res_ssvs_uniq$mcnv)
  print(table(res_ssvs_uniq$mcnv))
  
  print(paste0('Fraction mcnv-overlappers all: ', round((table1[2] / sum(table1)), 3)))
  print(paste0('Fraction mcnv-overlappers ssvs: ', round((table2[2] / sum(table2)), 3)))
  
  print(fisher.test(matrix(c(table1[1], table2[1], table1[2], table2[2]), nrow = 2), alternative='greater'))

}

overlap_with_cyto_bands <- function(res){

  script_link = '/Users/hoeps/PhD/projects/nahrcall/ntk-postprocess/scripts/bedtools_bands.sh'
  cyto_link = '/Users/hoeps/PhD/projects/nahrcall/analyses_paper/data/bands/cytoBand_simplified_renamed.bed'
  cyto_overlap_tmp_link = 'cyto_overlappers.bed'
  bedtools_link = '/usr/local/bin/bedtools'
  res_tmp_link = 'res_tmp.tmp'

  # Save res table, so bedtools can work with it.
  write.table(unique(res[,1:3]), res_tmp_link, sep='\t', col.names = F, row.names = F, quote = F)

  # Run bedtools to find intersections betwres_locus_meltn res and coredups
  bed_command = paste0(script_link ,' ', res_tmp_link, ' ', cyto_link,' > ', cyto_overlap_tmp_link)
  system(bed_command)

  # Load the results.
  res_uniq_cyto = read.table(cyto_overlap_tmp_link, sep='\t')
  colnames(res_uniq_cyto) = c('seqname','start','end', 'band')

  # Add this information to the original res directory
  res = dplyr::left_join(res, res_uniq_cyto, by=c('seqname', 'start','end'))
  res[is.na(res$band),'band'] = F

  return(res)
}

overlap_with_core_dups <- function(res){

  coredup_link = '/Users/hoeps/PhD/projects/huminvs/analyses_paper/data/genes/hg38/for_ntk/core-dup/cdups.bed'
  overlap_tmp_link = 'cd_overlappers.bed'
  bedtools_link = '/usr/local/bin/bedtools'
  res_tmp_link = 'res_tmp.tmp'

  # Save res table, so bedtools can work with it.
  write.table(res[,1:3], res_tmp_link, sep='\t', col.names = F, row.names = F, quote = F)

  # Run bedtools to find intersections betwres_locus_meltn res and coredups
  bed_command = paste0(bedtools_link ,' intersect -wa -a ', res_tmp_link, ' -b ', coredup_link, ' | sort | uniq > ', overlap_tmp_link)
  system(bed_command)

  # Load the results.
  cd_overlapping_segments = read.table(overlap_tmp_link, sep='\t')
  colnames(cd_overlapping_segments) = c('seqname','start','end')
  cd_overlapping_segments$core_dup = T

  # Add this information to the original res directory
  res = dplyr::left_join(res, cd_overlapping_segments, by=c('seqname', 'start','end'))
  res[is.na(res$core_dup),'core_dup'] = F

  return(res)
}

overlap_with_recurrent_invs <- function(res){

  rec_inv_link = '/Users/hoeps/PhD/projects/huminvs/analyses_paper/data/invs/recurrent/recurrent_invs.bed'
  overlap_tmp_link = 'inv_overlappers.bed'
  bedtools_link = '/usr/local/bin/bedtools'
  res_tmp_link = 'res_tmp.tmp'

  # Save res table, so bedtools can work with it.
  write.table(res[,1:3], res_tmp_link, sep='\t', col.names = F, row.names = F, quote = F)

  # Run bedtools to find intersections betwres_locus_meltn res and coredups
  bed_command = paste0(bedtools_link ,' intersect -wa -F 0.25 -a ', res_tmp_link, ' -b ', rec_inv_link, ' | sort | uniq > ', overlap_tmp_link)
  system(bed_command)

  # Load the results.
  inv_overlapping_segments = read.table(overlap_tmp_link, sep='\t')
  colnames(inv_overlapping_segments) = c('seqname','start','end')
  inv_overlapping_segments$rec_inv = T

  # Add this information to the original res directory
  res = dplyr::left_join(res, inv_overlapping_segments, by=c('seqname', 'start','end'))
  res[is.na(res$rec_inv),'rec_inv'] = F

  return(res)
}

overlap_with_mcnvs <- function(res){

  #mcnv_link = '/Users/hoeps/PhD/projects/nahrcall/analyses_paper/mcnv_enrichment/data/all_mcnvs_merged_plus50.bed'
  mcnv_link = '/Users/hoeps/PhD/projects/nahrcall/analyses_paper/mcnv_enrichment/data/all_mcnvs_merged_plus25pct.bed'
  overlap_tmp_link = 'mcnv_overlappers.bed'
  bedtools_link = '/usr/local/bin/bedtools'
  res_tmp_link = 'res_tmp.tmp'

  # Save res table, so bedtools can work with it.
  write.table(res[,1:3], res_tmp_link, sep='\t', col.names = F, row.names = F, quote = F)

  # Run bedtools to find intersections betwres_locus_meltn res and coredups
  bed_command = paste0(bedtools_link ,' intersect -wa -a ', res_tmp_link, ' -b ', mcnv_link, ' | sort | uniq > ', overlap_tmp_link)
  system(bed_command)

  # Load the results.
  mcnv_overlapping_segments = read.table(overlap_tmp_link, sep='\t')
  colnames(mcnv_overlapping_segments) = c('seqname','start','end')
  mcnv_overlapping_segments$mcnv = T

  # Add this information to the original res directory
  res = dplyr::left_join(res, mcnv_overlapping_segments, by=c('seqname', 'start','end'))
  res[is.na(res$mcnv),'mcnv'] = F

  return(res)
}

excl_effective_1_dup_del_pairs <- function(svdf){

  # Cut down to del/dup only
  svdf_deldup = svdf[svdf$sv %in% c('del','dup'),]

  svdf_deldup$len = svdf_deldup$end - svdf_deldup$start
  problematic = c()

  # Go through each pair (where coords partner2 > partner1).
  # For each, check if:
  # Case A: either their start or end matches, and their length differs by 1
  # Case B:
  for (n_partner1 in svdf_deldup$n){
    for (n_partner2 in svdf_deldup$n){
      svdf_pair = NULL
      if (n_partner2 > n_partner1){
        svdf_pair = svdf_deldup[svdf_deldup$n %in% c(n_partner1, n_partner2),]

        # Both have the same length?
        len_matches = abs(diff(svdf_pair$len)) <= 1

        # Both have same start coordinate?
        start_matches = diff(svdf_pair$start) == 0

        # Both have same end coordinate?
        end_matches = diff(svdf_pair$end) == 0

        # Both have same length but are offset by one (e.g. 1-3-dup  + 2-4-del).
        # This is not used as a criterion currently though.
        start_end_offset_one = (diff(svdf_pair$len) == 0) & (abs(diff(svdf_pair$start)) == 1)
        if (len_matches & (start_matches | end_matches )){#| start_end_offset_one)){
          problematic = c(problematic, svdf_pair$n)
        }
      }
    }
  }

  return(svdf[!(svdf$n %in% problematic),])
}



excl_shifter_inv_pairs <- function(svdf){
  svdf_inv = svdf[svdf$sv == 'inv',]

  i = svdf_inv$start
  j = svdf_inv$end
  xi = diff(svdf_inv$start)
  xj = diff(svdf_inv$end)
  #browser()
  # [10th Aug: was '<=1' originally, but I had chr7:144 which a inv-dup-inv where the two
  # Invs were on the same level.]
  problematic = svdf_inv[which((abs(xi + xj)) %in% c(1,2)),'n']

  return(svdf[!(svdf$n %in% c(problematic+1,problematic)),])
}

# Del dups that exist mh idk
excl_inverted_deldup <- function(svdf){

  inv_idx = which(svdf$sv == 'inv')
  del_idx = which(svdf$sv == 'del')

  if (length(del_idx) == 0 | (length(inv_idx) == 0)){
    return(svdf)
  }
  # If all dels come before invs, abort
  if (max(del_idx) < min(inv_idx)){
    return(svdf)
  }


  idx_swallowed = c()


  # For each inv, check if it overlaps with anything else.
  for (inv_n in inv_idx){
    inv_interval = data.frame(start = min(svdf[inv_n, c('start','end')]), end = max((svdf[inv_n, c('start','end')])))

    relevant_del_idx = del_idx[del_idx > inv_n]
    del_intervals = svdf[relevant_del_idx,]

    del_swallows = any((del_intervals$start <= inv_interval$start) & (del_intervals$end >= inv_interval$end))
    end_nested =   sum((svdf$end >= inv_interval$start) & (svdf$end <= inv_interval$end)) > 1


    if (del_swallows){
      idx_swallowed = c(idx_swallowed, inv_n)
    }
  }

  if (length(idx_swallowed ) > 0){
    return(svdf[-idx_swallowed,])
  } else {
    return(svdf)
  }
}

excl_swallowed_del <- function(svdf){

  del_idx = which(svdf$sv == 'del')

  if (length(del_idx) <= 1){
    return(svdf)
  }

  idx_swallowed = c()


  # For each del, check if it overlaps with another del.
  for (del_n in del_idx){
    del_interval = data.frame(start = min(svdf[del_n, c('start','end')]), end = max((svdf[del_n, c('start','end')])))

    relevant_del_idx = del_idx[del_idx > del_n]
    del_intervals = svdf[relevant_del_idx,]

    del_swallows = any((del_intervals$start <= del_interval$start) & (del_intervals$end >= del_interval$end))
    #end_nested =   sum((svdf$end >= inv_interval$start) & (svdf$end <= inv_interval$end)) > 1


    if (del_swallows){
      idx_swallowed = c(idx_swallowed, del_n)
    }
  }

  if (length(idx_swallowed ) > 0){
    return(svdf[-idx_swallowed,])
  } else {
    return(svdf)
  }
}


excl_swallowed_inv <- function(svdf){

  inv_idx = which(svdf$sv == 'inv')
  del_idx = which(svdf$sv == 'del')

  if (length(del_idx) == 0 | (length(inv_idx) == 0)){
    return(svdf)
  }
  # If all dels come before invs, abort
  if (max(del_idx) < min(inv_idx)){
    return(svdf)
  }


  idx_swallowed = c()


  # For each inv, check if it overlaps with anything else.
  for (inv_n in inv_idx){
    inv_interval = data.frame(start = min(svdf[inv_n, c('start','end')]), end = max((svdf[inv_n, c('start','end')])))

    relevant_del_idx = del_idx[del_idx > inv_n]
    del_intervals = svdf[relevant_del_idx,]

    del_swallows = any((del_intervals$start <= inv_interval$start) & (del_intervals$end >= inv_interval$end))
    end_nested =   sum((svdf$end >= inv_interval$start) & (svdf$end <= inv_interval$end)) > 1


    if (del_swallows){
      idx_swallowed = c(idx_swallowed, inv_n)
    }
  }

  if (length(idx_swallowed ) > 0){
    return(svdf[-idx_swallowed,])
  } else {
    return(svdf)
  }
}






filter_t <- function(t){
  #print(t)

  if (t == 'ref'){
    return('ref')
  }
  svdf = turn_mut_max_into_svdf(t, correction=T)

  # Filter 1: del can not swallow inv
  svdf = excl_swallowed_inv(svdf)

  # Filter 2: del can not swallow del
  svdf = excl_swallowed_del(svdf)
  # Filter 2: no 1-point dels
  svdf = svdf[!((svdf$sv == 'del') & (svdf$end - svdf$start == 1)),]
  svdf = svdf[!((svdf$sv == 'dup') & (svdf$end - svdf$start == 1)),]


  svdf = excl_effective_1_dup_del_pairs(svdf)

  # Filter 3: No i+x type double inversions
  svdf = excl_shifter_inv_pairs(svdf)

  # Transform back to string
  if(dim(svdf)[1] == 0){
    t_transform = 'ref'
  } else {
    t_transform = paste(paste(svdf$start, svdf$end, svdf$sv, sep='_'), collapse='+')
  }
  #print(t)
  #print(t_transform)
  #print('=========')
  return(t_transform)
}



# Annotate in the res file:
# mm0-mm4: this is mutation 0-4.
#
# In case of unresolved or contig break,
# m0 is Unresolved or Contig-break, and mm1-mm3 is NA or 'None'.
#
annotate_mut_better <- function(res_locus, solve_th = 98){

  res_locus$mut_max_cut = res_locus$mut_max
  res_locus = res_locus %>% separate(mut_max_cut, c("mm1", "mm2", "mm3"), sep="\\+")
  res_locus$mm0 = 'ref'
  res_locus[res_locus$mm1 == 'ref','mm1'] = NA


  # If size_mut > 1 and there is no inv overlap, category 'Non-nested'
  #res_locus$overlapthing = res_locus$inv_overlap + res_locus$any_overlap

  # autodetect overlap mode
  res_locus[(res_locus$size_mut > 1) & (res_locus$overlap == F),'mm0'] = 'Non-nested'
  res_locus[res_locus$res_max < solve_th,'mm0'] = 'Unresolved'
  if (any(res_locus$exceeds_y==T)){
    res_locus[res_locus$exceeds_y==T,'mm0'] = 'Contig-break'
  }
  res_locus[res_locus$res_max < solve_th,c('mm1', 'mm2', 'mm3')] = 'None'
  res_locus[res_locus$exceeds_y==T,c('mm1', 'mm2', 'mm3')] = 'None'
  res_locus[(res_locus$size_mut > 1) & (res_locus$overlap == F),c('mm1', 'mm2', 'mm3')] = 'None'


  # [W, 15th Aug 2022]
  # Here is some filtering stuff going on, BUT we filter already with filter_t
  # so I think this is not appropriate anymore.

  #####
  # exists = any(na.omit(res_locus$mm3 == res_locus$mm1))
  #
  # if (exists){
  #   res_locus[res_locus$mm3 == res_locus$mm1,]$mm3 = NA
  # }
  # res_locus[res_locus$mm2 == res_locus$mm1,]$mm2 = NA
  # if (any(res_locus$mm1 == res_locus$mm0)){
  #   res_locus[res_locus$mm1 == res_locus$mm0,]$mm1 = NA
  # }
  #####

  res_locus$mm1 = gsub("^.*_.*_","",res_locus$mm1)
  res_locus$mm2 = gsub("^.*_.*_","",res_locus$mm2)
  res_locus$mm3 = gsub("^.*_.*_","",res_locus$mm3)


  res_locus[res_locus$result=='Not explained',c('mm1','mm2','mm3')] = NA

  if (dim(res_locus[res_locus$result=='Not explained' & res_locus$mm0=='ref',])[1] > 0){
    res_locus[res_locus$result=='Not explained' & res_locus$mm0=='ref','mm0'] = 'Unresolved'
  }

  if (dim(res_locus[res_locus$result=='Not explained' & res_locus$mm0=='None',])[1] > 0){
    res_locus[res_locus$result=='Not explained' & res_locus$mm0=='None','mm0'] = 'Unresolved'
  }

  return(res_locus)
}

turn_mut_max_into_svdf <- function(t, correction=T){
  if (t == 'ref'){
    return(F)
  }
  
  
  strsplit(t, split='n+')
  sv_list = strsplit(t, split="\\+")[[1]]
  svdf_clump = data.frame(sv_list)
  svdf = data.frame(do.call('rbind', strsplit(as.character(svdf_clump$sv_list),'_',fixed=TRUE)))
  colnames(svdf) = c('start', 'end', 'sv')

  svdf$start = as.numeric(svdf$start)
  svdf$end = as.numeric(svdf$end)
  svdf$n = 1: dim(svdf)[1]


  if (correction == F){
    return(svdf)
  }




  # Iteratve over each sv
  for (i in 1:dim(svdf)[1]){
    
    
    # if it's an inv, complicated changes have to happen.
    if (svdf[i,'sv'] == 'inv'){

      # We work only on things that are inside of the inversion. Stuff that partially
      # overlaps it are simply too complicated to deal with.

      #end breakpoints get mirrored
      svdf[(svdf$start >= svdf[i,'start']) & (svdf$end <= svdf[i,'end']) & (svdf$n > i),'newend'] =
        svdf[i,'end'] + (svdf[(svdf$start >= svdf[i,'start']) & (svdf$end <= svdf[i,'end']) & (svdf$n > i),'start'] - svdf[i,'start'])

      #start breakpoints get mirrored
      svdf[(svdf$start >= svdf[i,'start']) & (svdf$end <= svdf[i,'end']) & (svdf$n > i),'newstart'] =
        svdf[i,'start'] + (svdf[i,'end'] - svdf[(svdf$start >= svdf[i,'start']) & (svdf$end <= svdf[i,'end']) & (svdf$n > i),'end'] )
    }

    svdf[!is.na(svdf$newend),'end'] = svdf[!is.na(svdf$newend),'newend']
    svdf[!is.na(svdf$newstart),'start'] = svdf[!is.na(svdf$newstart),'newstart']
    svdf[,c('newend', 'newstart')] = NULL
    
    # If it's del/dup, change the del coordinates.
    if (svdf[i,'sv'] == 'del'){

      #start breakpoints after the beginning of del get a plus
      svdf[(svdf$start > svdf[i,'start']) & (svdf$n > i),'start'] =
        svdf[(svdf$start > svdf[i,'start']) & (svdf$n > i),'start'] + (svdf[i,'end'] - svdf[i,'start'])

      #end breakpoints after the beginning of del get a plus
      svdf[(svdf$end > svdf[i,'start']) & (svdf$n > i),'end'] =
        svdf[(svdf$end > svdf[i,'start']) & (svdf$n > i),'end'] + (svdf[i,'end'] - svdf[i,'start'])

    } else if (svdf[i,'sv'] == 'dup'){

      # start breakpoints after the dup get a minus
      svdf[(svdf$start > svdf[i,'end']) & (svdf$n > i),'start'] =
        svdf[(svdf$start > svdf[i,'end']) & (svdf$n > i),'start'] - (svdf[i,'end'] - svdf[i,'start'])

      # end breakpoints after thhe dup get a minus
      svdf[(svdf$end > svdf[i,'end']) & (svdf$n > i),'end'] =
        svdf[(svdf$end > svdf[i,'end']) & (svdf$n > i),'end'] - (svdf[i,'end'] - svdf[i,'start'])
    }
  }
  return (svdf)
}


determine_inv_overlapper <- function(t, overlap_mode){

  if (t == 'ref'){
    return(F)
  }

  svdf = turn_mut_max_into_svdf(t, correction=T)

  # No inv, no inv overlapper...
  if ((overlap_mode=='inv_overlap') & (!('inv' %in% svdf$sv))){
    return(F)
  }

  if (overlap_mode == 'inv_overlap'){
    dfi = svdf[svdf$sv == 'inv',]
    any_overlap = F

    # For each inv, check if it overlaps with anything else.
    for (inv_n in 1:dim(dfi)[1]){
      inv_interval = data.frame(start = min(dfi[inv_n, c('start','end')]), end = max((dfi[inv_n, c('start','end')])))
      start_nested = sum((svdf$start >= inv_interval$start) & (svdf$start <= inv_interval$end)) > 1
      end_nested =   sum((svdf$end >= inv_interval$start) & (svdf$end <= inv_interval$end)) > 1
      inv_nested = sum((svdf$end >= inv_interval$end) & (svdf$start <= inv_interval$start)) > 1

      if (start_nested | end_nested | inv_nested){
        any_overlap = T
      }
    }
  } else if (overlap_mode == 'any_overlap'){

    any_overlap = F
    for (sv_n in 1:dim(svdf)[1]){
      interval = data.frame(start = min(svdf[sv_n, c('start','end')]), end = max((svdf[sv_n, c('start','end')])))
      start_nested = sum((svdf$start >= interval$start) & (svdf$start <= interval$end)) > 1
      end_nested =   sum((svdf$end >= interval$start) & (svdf$end <= interval$end)) > 1


      if (start_nested | end_nested){
        any_overlap = T
      }
    }
  }
  # print(svdfo)
  # print(svdf)
  # print(any_overlap)
  return(any_overlap)
}


excuse_missing_samples <- function(res){

  r = unique(res)
  n_samples = max(table(r$start))
  for (start in unique(r$start)){
    #print(start)
    q = r[r$start==start,]
    q_model = q[1,]
    missing_samples = unique(r$sample)[!unique(r$sample) %in% q$sample]

    for (m in missing_samples){
      new_col = q_model[1,]
      new_col$sample = m
      new_col[,c('res_ref', 'res_max', 'mut_max', 'mut_simulated',
                 'mut_tested', 'search_depth', 'grid_compression', 'flip_unsure'
      )] = c(0,0,'ref', 0,
             0,0,0,F)
      q_model = rbind(q_model, new_col)
    }
    r = rbind(r, q_model)
  }

  r$res_ref = as.numeric(r$res_ref)
  r$res_max = as.numeric(r$res_max)
  r$mut_simulated = as.numeric(r$mut_simulated)
  r$mut_tested = as.numeric(r$mut_tested)
  r$search_depth = as.numeric(r$search_depth)
  r$grid_compression = as.numeric(r$grid_compression)

  return(unique(r))
}



find_reslink <- function(twenty=F){

  res_link = 'data/above50_8th_apr.tsv'
  res_link = '~/PhD/projects/nahrcall/all_samps_may12.tsv'
  res_link = '~/Desktop/pics_fullrun_jul11/all.tsv'
  res_link = '/Users/hoeps/PhD/projects/nahrcall/ntk-postprocess/res_15_jul/all.tsv'
  res_link = '/Users/hoeps/PhD/projects/nahrcall/ntk-postprocess/res_15_jul/all_plusnachzug.tsv'
  #res_link = '/Users/hoeps/PhD/projects/nahrcall/analyses_paper/cao_regions/res/all.tsv'
  #res_link = '/Users/hoeps/PhD/projects/nahrcall/analyses_paper/17-interesting/calls/results_23.tsv'
  #res_link = '/Users/hoeps/PhD/projects/nahrcall/analyses_paper/17-interesting/calls/calls_errorfixing_22_jul/all_uniq.tsv'
  #res_link = '/Users/hoeps/PhD/projects/nahrcall/analyses_paper/17-interesting/calls/calls_errorfixing_22_jul/all_uniq_7_exchanged_uniq.tsv'
  #res_link = '/Users/hoeps/PhD/projects/nahrcall/analyses_paper/17-interesting/calls/calls_25th.tsv'

  res_link = '/Users/hoeps/PhD/projects/nahrcall/analyses_paper/inv_full_run_1st_aug/calls/all_u.tsv'
  res_link = '/Users/hoeps/PhD/projects/nahrcall/analyses_paper/17-interesting/calls/calls_8_aug.tsv'
  res_link = '/Users/hoeps/PhD/projects/nahrcall/analyses_paper/17-interesting/calls/calls_aug_9/all_plus7.tsv'
  res_link = '/Users/hoeps/PhD/projects/nahrcall/analyses_paper/17-interesting/calls/calls_aug_12_98percent/all.tsv'
  res_link = '/Users/hoeps/PhD/projects/nahrcall/analyses_paper/17-interesting/calls/calls_15aug/all.tsv'
  res_link = '/Users/hoeps/PhD/projects/nahrcall/analyses_paper/17-interesting/calls/calls_15aug/all_plus_chm13.tsv'
  res_link = '/Users/hoeps/PhD/projects/nahrcall/analyses_paper/17-interesting/calls/calls_15aug/integrate_chm13/all_plus_chm13.tsv'
  res_link = '/Users/hoeps/PhD/projects/nahrcall/analyses_paper/17-interesting/calls/calls_15aug/integrate_chm13/integrate_new_chr15/all_all_chm13_aug23.tsv'
  res_link = '/Users/hoeps/PhD/projects/nahrcall/analyses_paper/17-interesting/calls/calls_15aug/integrate_chm13/integrate_new_chr15/exchange_chr22/all3_aug25.tsv'
  
  #res_link = '/Users/hoeps/PhD/projects/nahrcall/analyses_paper_2/first_run/all.tsv'
  res_link = '/Users/hoeps/PhD/projects/nahrcall/analyses_paper_2/second_run/all.tsv'
  
  #res_link = '/Users/hoeps/PhD/projects/nahrcall/analyses_paper/figures/Fig3/sotos/ship/calls_origonly.tsv'

  if (twenty){
    res_link = '/Users/hoeps/PhD/projects/nahrcall/analyses_paper/17-interesting/calls/calls_morres_locus_meltrrorfix_27_jul/all_u.tsv'
  }
  return(res_link)
}

make_dumres_locus_m_onlyell_plot <- function(res1samp){


  res1samp$invname = paste(res1samp$seqname, res1samp$start, res1samp$end, sep = '-')

  res1samp$simple_sample =
    sub("_NA", "",
        paste0(
          sub("_[hg].*", "", res1samp$sample),
          '_',
          str_extract(res1samp$sample, "h[12]")
        )
    )


  # Sort so the apes are in the bottom
  apenames = c('panPan3', 'panTro6','gorGor6','ponAbe3','rheMac10')
  nonapes = res1samp$simple_sample[!res1samp$simple_sample %in% apenames]
  better_order = c(sort(nonapes), apenames)

  res1samp = res1samp %>%
    slice(match(better_order, simple_sample))
  res1samp$n_line = dim(res1samp)[1]:1


  p = ggplot(data=res1samp) + geom_segment(aes(x=res_ref, xend=res_max, y=n_line, yend=n_line, color=exceeds_y)) +
    geom_point(aes(x=res_max, y=n_line), color='red') +
    geom_point(aes(x=res_ref, y=n_line), color='black') +
    geom_point(data=res1samp[res1samp$exceeds_y==T,], aes(x=res_ref, y=n_line), color='red') +
    geom_text(aes(x=105, y=n_line, label= paste0(simple_sample, ' ', mut_max)), size=3, hjust=0) +
    theme_bw() +
    scale_x_continuous(breaks=seq(0,100,10), limits=c(0,200)) +
    labs(x = 'Sequence concordance with hg38',
         title = paste0(res1samp[1,'invname']))

  p
  return(p)
}

res_find_sSVs <- function(res_f, overlap_mode, ape_samples, length_limit_bp = 20000){

  res_f_plus = res_f

  # Filter
  res_f_plus_filt = res_f_plus[(res_f_plus$exceeds_y == F) &
                                 !(res_f_plus$sample %in% ape_samples) &
                                 ((res_f_plus$end - res_f_plus$start) > length_limit_bp),]

  # What was cured by >1 mutation?
  res_cure = res_f_plus_filt[(res_f_plus_filt$cured==T) & (res_f_plus_filt$size_mut > 1),]

  if (overlap_mode == 'inv_overlap'){
    # ... and has an inv?
    res_cure_interest = res_cure[(res_cure$onlydels + res_cure$onlydups == 0) & (res_cure$inv_involved == T),]
  } else {
    res_cure_interest = res_cure
  }
  # ... and overlaps with sth
  res_cure_interest$overlap =  unlist(lapply(res_cure_interest$mut_max, determine_inv_overlapper, overlap_mode = overlap_mode))
  res_cure_interest = res_cure_interest[res_cure_interest$overlap==T,]


  return(res_cure_interest)
}

add_sizemut_invinvolved_etal_to_res <- function(res_f, cure_threshold_pct = 98){
  
  res_f$mut_max = unlist(lapply(res_f$mut_max, filter_t))
  
  res_f = excuse_missing_samples(res_f)

  # Info #1
  res_f$size_mut = lengths(regmatches(res_f$mut_max, gregexpr("\\+", res_f$mut_max))) + 1

  # Info #2
  res_f$inv_involved = lengths(regmatches(res_f$mut_max, gregexpr("inv", res_f$mut_max))) > 0

  # Info #3
  res_f[res_f$mut_max == 'ref',]$size_mut = 0

  # Info #4
  res_f$cured = res_f$res_max >= cure_threshold_pct
  
  # Info #5
  res_f$n_inv_involved = str_count(res_f$mut_max, 'inv')
  res_f$n_del_involved = str_count(res_f$mut_max, 'del')
  res_f$n_dup_involved = str_count(res_f$mut_max, 'dup')

  # Info #6
  res_f$inv_involved = res_f$n_inv_involved > 0

  # Info #7
  res_f$onlydels = res_f$n_del_involved == res_f$mut_max
  res_f$onlydups = res_f$n_dup_involved == res_f$mut_max

  # Info #8
  res_f$res_fult = 'NA'
  res_f[res_f$res_max < cure_threshold_pct, 'result'] = 'Not explained'
  res_f[res_f$cured==T, 'result'] = res_f[res_f$cured==T, 'size_mut']
  #res_f[res_f$res_ref >= cure_threshold_pct, 'result'] = 'Ref'
  res_f[res_f$result %in% c('4','5','6','7','8','9','10'), 'result'] = "4+"
  res_f[res_f$exceeds_y==T,'result'] = 'Contig-break'

  return(res_f)
}

plot_one_locus_v1 <- function(res, start){

  res_one_locus = res[res$start == start,]

  res_one_locus = within(res_one_locus, mm<-data.frame(do.call('rbind', strsplit(as.character(res_one_locus$mut_max), '+', fixed=TRUE))))
  res_one_locus$mm0 = 'ref'
  res_one_locus$mm1 = res_one_locus$mm$X1
  res_one_locus$mm2 = res_one_locus$mm$X2
  res_one_locus$mm3 = res_one_locus$mm$X3

  # res_one_locus[res_one_locus$mm3 == res_one_locus$mm2,]$mm3 = NA
  # res_one_locus[res_one_locus$mm2 == res_one_locus$mm1,]$mm2 = NA
  # res_one_locus[res_one_locus$mm1 == res_one_locus$mm0,]$mm1 = NA

  res_one_locus$mm1 = gsub("^.*_.*_","",res_one_locus$mm1)
  res_one_locus$mm2 = gsub("^.*_.*_","",res_one_locus$mm2)
  res_one_locus$mm3 = gsub("^.*_.*_","",res_one_locus$mm3)

  res_mm_only = res_one_locus[,c('sample', 'mm0','mm1', 'mm2', 'mm3')]
  res_mm_molten = reshape2::melt(res_mm_only,  id.vars = 'sample')

  dd = res_mm_molten[order(res_mm_molten$value),]
  p = ggplot(dd) + geom_tile(aes(y=sample, x=variable, fill=value))

  print(p)
}


plot_one_locus_v2 <- function(res_f, anc, start_coord, solve_th = 98){


  res_f = annotate_mut_better(res_f)
  # choose sample
  res_locus = res_f[res_f$start == start_coord,]
  #res_locus = annotate_mut_better(res_locus, solve_th = solve_th)

  # Process
  res_locus_m_only = res_locus[,c('sample','mm0', 'mm1', 'mm2', 'mm3')]

  # Add ancestry information
  res_locus_m_only$simplesample = sub("\\_.*","", sub("\\..*", "", res_locus_m_only$sample))
  res_locus_m_only = left_join(res_locus_m_only, anc, by='simplesample')
  res_locus_m_only =   within(res_locus_m_only,    mm0 <-    factor(mm0, levels=c('ref','Contig-break','Unresolved')))
  
  res_locus_m_only = res_locus_m_only[with(res_locus_m_only, order(mm0, mm1,mm2, mm3, ANC)),]
  
  # Determine x order
  y_order = data.frame(sample = res_locus_m_only$sample, n = 1:length(res_locus_m_only$sample))

  res_locus_melt = reshape2::melt(res_locus_m_only,  id.vars = c('sample','simplesample'))
  res_locus_melt = res_locus_melt[order(res_locus_melt$value),]
  res_locus_melt = left_join(res_locus_melt, y_order, by='sample')


  res_locus_melt$width = 1
  widthval = 0.5
  res_locus_melt[res_locus_melt$variable=='mm0', ]$width = widthval
  res_locus_melt[res_locus_melt$variable=='mm1', ]$width = widthval
  res_locus_melt[res_locus_melt$variable=='mm2', ]$width = widthval
  res_locus_melt[res_locus_melt$variable=='mm3', ]$width = widthval
  
  res_locus_melt[is.na(res_locus_melt$value),'value'] = 'NA'
  

  color_vals = list(
    'AMR' = '#671128',
    'AFR' = '#F7D462',
    'EAS' = '#7B8028',
    'EUR' = '#3E89A7',
    'SAS' = '#8663A7',
    'hg38' = '#E86BA5',
    'Unresolved' = '#D3D3D3',
    'Contig-break' = '#BEBDBD',
    'ref' = '#514B47',
    'inv' = '#4D9E47',
    'del' = '#D14D2A',
    'dup' = '#529FD5',
     'NA' = '#FFFFFF'
  )
  

  color_vals = color_vals[names(color_vals) %in% unique(res_locus_melt$value)]
  
  res_locus_melt$variable <- with(res_locus_melt,factor(variable,levels = c('ANC','mm0','mm1','mm2','mm3')))
  res_locus_melt$plotcolor = as.character(color_vals[res_locus_melt$value])
  res_locus_melt =   within(res_locus_melt,    value <-    factor(value, levels=names(color_vals)))
  res_locus_melt[is.na(res_locus_melt$value),'value'] = 'NA'
  
  outplot = ggplot(res_locus_melt) + 
    geom_tile(aes(x=variable, y=n, fill=value, width=width), color='black', height=1) +
    scale_fill_manual(values=as.character(color_vals)) +
    scale_y_discrete(labels = sample, breaks=1:length(res_locus_melt$sample)) +
    scale_x_discrete(labels = c('Ancestry','Ref','Mut #1','Mut #2','Mut #3')) +
    labs(title=start_coord, x='')

  print(outplot)

  return(list(outplot,res_locus))
}


plot_overview_n_mut <- function(res_f, sort_by_ref=F){

  res_f = annotate_mut_better(res_f)
  res_f[res_f$mm0 == 'Non-nested', 'result'] = 'Non-Nested'
  res_sum = (res_f %>% group_by(start) %>% mutate(failed_call = sum(result=='Not explained') + sum(result=='Contig-break')) %>% slice(1))
  #res_sum = (res_f %>% group_by(start) %>% mutate(failed_call = sum(result=='Contig-break')) %>% slice(1))
  res_sum2 = (res_f %>% group_by(start) %>% mutate(failed_call_ref = sum(result=='Not explained') + sum(result=='Contig-break') + sum(result=='0')) %>% slice(1))

  if (sort_by_ref){
    y_order = res_sum2[order(res_sum2$failed_call_ref),c('seqname','start', 'end')]
  } else {
    y_order = res_sum[order(res_sum$failed_call, res_sum$start),c('seqname','start', 'end')]
  }

  #browser()
  start_order = y_order$start
  start_order = as.character(as.data.frame(start_order)$start)
  fillorder = rev(c( 'Contig-break','Not explained','Non-Nested','0','1','2','3', '4+'))
  #cbPalette <- rev(c("#b342f5", "#d14f56", "#6c79b8", "#bfb743", "#1dc200", "#949494", 'black'))

  color_seven = (length(unique(res_f$result)) == 7)

  if (color_seven){
    cbPalette <- rev(c('grey', 'lightgrey','lightblue','green', 'yellow','orange','red'))

  } else {
    cbPalette <- rev(c('grey', 'lightgrey','lightblue','green', 'yellow','orange','red', 'black'))
}
  res_f =   within(res_f,    resi <-          factor(result, levels=fillorder))
  res_f =   within(res_f,    start_fact <-    factor(start, levels=rev(start_order)))




  chr_band_dict = unique(res_f[,c('seqname','start','end','band', 'width_orig')])
  chr_band_dict$displayed_name = paste0(chr_band_dict$band, ' (', round(chr_band_dict$width_orig / 1000, 0), 'kbp)')
  row.names(chr_band_dict) = chr_band_dict$start#paste(chr_band_dict$seqname, chr_band_dict$start, chr_band_dict$end, sep='-')
  sorted_display_names = chr_band_dict[levels(res_f$start_fact),'displayed_name']


  p = ggplot(res_f) + geom_bar(aes(y=start_fact, fill=resi)) +
    scale_fill_manual(values=cbPalette, name='# Sequential Mutations') +
    theme_bw() +

    scale_y_discrete(labels = sorted_display_names) +
    scale_x_continuous(breaks = c(seq(0,length(unique(res_f$sample)), 20), length(unique(res_f$sample)), 20)) +
    labs(x='Number of samples', y='Genomic locus')

  print(p)

  return(list(p, start_order))
}


table_sort_and_prep <- function(table, sv_orders_link){
  svo = read.table(sv_orders_link, header=F, sep='\t')
  colnames(svo) = 'mut'

  svo$svo = (dim(svo)[1] + 1) - as.numeric(row.names(svo))

  row.names(table) = str_replace(
    str_replace(
      (str_replace_all(row.names(table), " NA", "")),
      'ref ', ''
    ),
    ' None', '')
  row.names(table)[row.names(table) == 'Non-nested None None'] = 'Non-nested'

  table$mut = row.names(table)

  table2 = left_join(table, svo, by='mut')


  row.names(table2) = table2$mut
  table2$mut = NULL

  table2['Resolved',] = colSums(table2) - colSums(table2[c('Unresolved', 'Contig-break', 'napes', 'Core_duplicon'),])
  table2['Resolved','svo'] = 1
  table2 = table2[order(table2$svo, decreasing = F),]

  table2$svo = NULL

  return(table2)
}

second.word <- function(my.string){
  unlist(strsplit(my.string, "-"))[2]
}

make_inferno_heatmap <- function(res_plus, solve_th, process_all = F){

  sv_orders_link = '/Users/hoeps/PhD/projects/nahrcall/analyses_paper/17-interesting/plot/sv_orders2.txt'
  # Makes resp2, a subset of 'res' which contains only locations in which ssvs have bres_locus_meltn sres_locus_meltn.
  # if (process_all == F){
  #   resp2 = res_plus[res_plus$start %in% unique(ssvs$start),]
  # } else if (process_all == T){
  # }
  resp2 = res_plus
  
  # Add a whole lot of annotation to resp2.
  resp2 = annotate_mut_better(resp2, solve_th = solve_th)
  resp2[resp2$mm0=='Contig-break',c('mm1','mm2','mm3')] = NA
  resp2$mut_plot = (paste(resp2$mm0, resp2$mm1, resp2$mm2, resp2$mm3))
  resp2$region= paste(resp2$seqname, resp2$start, resp2$end, sep='-')

  if (process_all){
    return(list(resp2, resp2))
  }

  # Make a table.
  table = as.data.frame(t(as.data.frame.matrix(table(resp2$region, resp2$mut_plot))))

  # Add core-dups and apes post-hoc
  for (add_column in c('n_ape_resolved','core_dup','rec_inv', 'mcnv')){
    table = add_datacol_to_table(table, resp2, add_column)
  }
  # add_datacol_to_table(table, resp2, 'n_ape_resolved')
  #
  # table = add_apes_to_table(table, resp2, 'n_ape_resolved')
  # table = add_core_dups_to_table(table, resp2, 'core_dup')
  # table = add_invs_to_t dable(table, resp2, 'rec_inv')

  # Sort and prep.
  table2 = table_sort_and_prep(table, sv_orders_link)
  table_tmp = table2

  table_tmp['failed_call',] = table_tmp['Unresolved',] + table_tmp['Contig-break',]
  table_tmp['start',] = as.numeric(sapply(colnames(table_tmp), second.word))
  # table2 = table2[,order(as.numeric(table2['Unresolved',] + table2['Contig-break',]))]
  # table2 = table2[,order(as.numeric(table2['Contig-break',]))]

  table2 = table2[,order(as.numeric(table_tmp['failed_call',]), as.numeric(table_tmp['start',]), decreasing = F)]

  # Replace chr-start-end with band + length. Ein kleiner Exkurs...
  chr_band_dict = unique(resp2[,c('seqname','start','end','band', 'width_orig')])
  chr_band_dict$displayed_name = paste0(chr_band_dict$band, ' (', round(chr_band_dict$width_orig / 1000, 0), 'kbp)')
  row.names(chr_band_dict) = paste(chr_band_dict$seqname, chr_band_dict$start, chr_band_dict$end, sep='-')
  colnames(table2) = chr_band_dict[colnames(table2),'displayed_name']
  #Exkurs over

  #colnames(table2) = unique(paste0(resp2$band, ' (', round(resp2$width_orig / 1000, 0), 'kbp)'))

  mat_breaks = c(0,0.25,0.5,0.75,1,1.5,2,3,4,5,10,20,40, 58)

  t3 = t(table2)
  t4 = cbind(t3[,6:8], t3[,1:5], t3[,9:ncol(t3)])
  t4[,1:3] = t4[,1:3] * -1

  p = pheatmap((t4),
               #color  = c("#FFFFFF", inferno(length(mat_breaks)-1, direction=1)),
               color  = c("#4dd9ff","#000000", inferno(length(mat_breaks), direction=1)),

               cluster_cols = F,
               cluster_rows = F,
               border_color = 'black',
               breaks  = c(-1,-0.1, mat_breaks),
               cutree_cols = 3,
               gaps_col = c(3,3,8,19))



  print(p)
  return(list(p,resp2))
}




add_n_ape_resolved <- function(res_f, ape_samples, res_th = 98){

  res_ape = res_f[res_f$sample %in% ape_samples,]
  res_ape$ape_resolve = F
  res_ape[res_ape$res_max > res_th,]$ape_resolve = T
  n_resolve = res_ape %>% group_by(start) %>% mutate(sum(ape_resolve==T)) %>% slice(1)
  n_resolve_sliced = res_ape %>% group_by(start) %>% mutate(n_ape_resolved = sum(ape_resolve==T)) %>% slice(1)

  n_resolve_sliced_condensed = n_resolve_sliced[,c('start','n_ape_resolved')]
  res_f_apes = left_join(res_f, n_resolve_sliced_condensed, by='start')

  return(res_f_apes)
}

add_datacol_to_table <- function(table, resp2, colname){

  if (colname == 'n_ape_resolved'){
    added_df = data.frame(region=resp2$region, napes = resp2[[colname]]) %>% group_by(region) %>% slice(1)
  } else if (colname == 'core_dup'){
    added_df = data.frame(region=resp2$region, Core_duplicon = as.numeric(resp2[[colname]])) %>% group_by(region) %>% slice(1)
  } else if (colname == 'rec_inv'){
    added_df = data.frame(region=resp2$region, Recurrent_Inv = as.numeric(resp2[[colname]])) %>% group_by(region) %>% slice(1)
  } else if (colname == 'mcnv'){
    added_df = data.frame(region=resp2$region, mCNV = as.numeric(resp2[[colname]])) %>% group_by(region) %>% slice(1)
  }

  t_table2 = as.data.frame(t(table))
  t_table2$region = row.names(t_table2)
  t_table3 = left_join(t_table2, added_df, by='region')
  row.names(t_table3) = t_table3$region
  t_table3$region = NULL
  t_table_out = t(t_table3)

  return(as.data.frame(t_table_out))
}

# add_core_dups_to_table <- function(table, resp2, colname){
#   added_df = data.frame(region=resp2$region, Core_duplicon = as.numeric(resp2[[colname]])) %>% group_by(region) %>% slice(1)
#
#   t_table2 = as.data.frame(t(table))
#   t_table2$region = row.names(t_table2)
#   t_table3 = left_join(t_table2, added_df, by='region')
#   row.names(t_table3) = t_table3$region
#   t_table3$region = NULL
#   t_table_out = t(t_table3)
#
#   return(data.frame(t_table_out))
# }

evaluate <- function(res, res_test){
  curated_intervals_link = '/Users/hoeps/PhD/projects/nahrcall/analyses_paper/17-interesting/23_regions.bed'
  ci = read.table(curated_intervals_link, sep='-')
  colnames(ci) = c('seqnames','start','end')
  ci = ci[ci$start %in% res$start,]

  frac_retained = sum(ci$start %in% res_test$start) / length(ci$start)

  retained = ci[ci$start %in% res_test$start,]
  not_retained = ci[!(ci$start %in% res_test$start),]

  print(paste0('Refined/Unrefined fraction: ', round(length(unique(ssvs$start))/length(ci$start),3)))
  print(paste0('Retained fraction: ', round(frac_retained,3)))
  print('Missing:')
  print(not_retained)
}

res_mut_into_sep_columns <- function(res){
  res_f = within(res_f, mm<-data.frame(do.call('rbind', strsplit(as.character(res_f$mut_max), '+', fixed=TRUE))))

  res_f$mm0 = 'ref'
  res_f$mm1 = res_f$mm$X1
  res_f$mm2 = res_f$mm$X2
  res_f$mm3 = res_f$mm$X3

  res_f[res_f$mm3 == res_f$mm1,]$mm3 = NA
  res_f[res_f$mm2 == res_f$mm1,]$mm2 = NA
  res_f[res_f$mm1 == res_f$mm0,]$mm1 = NA

  res_f = within(res_f, m1<-data.frame(do.call('rbind', strsplit(as.character(res_f$mm1), '_', fixed=TRUE))))
  res_f = within(res_f, m2<-data.frame(do.call('rbind', strsplit(as.character(res_f$mm2), '_', fixed=TRUE))))
  res_f = within(res_f, m3<-data.frame(do.call('rbind', strsplit(as.character(res_f$mm3), '_', fixed=TRUE))))

  res_f$m1start = res_f$m1$X1
  res_f$m1end = res_f$m1$X2
  res_f$m1name = res_f$m1$X3

  res_f$m2start = res_f$m2$X1
  res_f$m2end = res_f$m2$X2
  res_f$m2name = res_f$m2$X3

  res_f$m3start = res_f$m3$X1
  res_f$m3end = res_f$m3$X2
  res_f$m3name = res_f$m3$X3

  res_f[,c('mm1','mm2','mm3', 'mm', 'm1','m2','m3')] = NULL

  return(res_f)
}

copy_ssvs_pdf_to_local <- function(ssvs_full, link_to_remote_res, link_to_local_target, out_script, execute=F){
  
  ssvs_full_process = ssvs_full
  ssvs_full_process$id <- paste(ssvs_full_process$seqname, ssvs_full_process$start, ssvs_full_process$end, sep="-")
  
  
  id_cmds = paste0('rsync -R ', 
                   link_to_remote_res, '/', unique(ssvs_full_process$id), '/*/*/*.pdf ', 
                   link_to_remote_res, '/', unique(ssvs_full_process$id), '/*/*/*/*.pdf ', 
                   link_to_local_target)
  if (execute == T){
    print('Copying ssv pdfs over...')
    cp_count = 0
    for (cmd in id_cmds){
      cp_count = cp_count + 1
      print(paste0('Copying folder ', cp_count, ' out of ', length(id_cmds)))
      system(cmd)
    }
  } else {
    print('Not executing the cp. Writing commands instead.')
    write.table(id_cmds, file=out_script, row.names=F, col.names=F, quote=F)
    print(paste0('Bashfile for copying written to: ', out_script))
  }
  
}


return_supp_ufo_higher <- function(res, params){
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
  heat = Heatmap(rescast2[,sort_order], col = colors, cluster_columns = F, cluster_rows = F,
                 column_split = rep(1:5, as.numeric(table(res_slice$category))))
  return(heat)
}

return_supp_ufo_lower <- function(res, params){
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
  pheat = pheatmap(res_plot,
                   cluster_cols = F,
                   cluster_rows = F,
                   gaps_col = white_breaks,
                   color = c('#2d2e30', hcl.colors(10, "BluYl")))
  
  return(list(res_slice, pheat))
}
# Make an overview chart

return_fig2a_bars <- function(res_slice){
  cat = as.data.frame(table(res_slice$category))
  cat$Var1 = as.character(cat$Var1)
  cat[cat$Var1 == 'No Contig ASM', 'Var1'] = 'No contiguous asm'
  cat$Var1 = factor(cat$Var1, levels=c('No contiguous asm', 'Unexplained SVs', 'Ref','Simple SVs','sSVs'))
  
  cols = c('#2d2e30', '#585b61', '#466299', '#8eba7f', '#edb03e', '#ed583e')
  
  
  p = ggplot(data = cat) + geom_bar(aes(x=Var1, y=Freq, fill=Var1), stat = 'identity') + theme_bw() + 
    labs(x='Locus group', y='# Loci') + 
    scale_fill_manual(values = cols)
  
  return(p)
  
  #ggsave(p, file = '../plots/FigOverviewbars.pdf', device='pdf', width = 14, height = 10 , units='cm')
}

