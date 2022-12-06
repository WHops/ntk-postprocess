# Laboratory


t1 = '16_33_inv+33_41_del+4_10_del'
t2 = '16_44_inv+12_28_del'
t3 = '14_18_dup+3_18_inv+16_19_del'
t4 = '1_11_dup+5_8_del+12_15_del+13_18_inv+19_21_del'
t5 = '4_6_dup+7_12_inv'
t6 = '2_15_dup+7_16_inv'
t7 ='9_19_inv+19_24_dup+19_23_del'
t = t7

determine_inv_overlapper(t6)

determine_inv_overlapper <- function(t){
  strsplit(t, split='+')
  sv_list = strsplit(t, split="\\+")[[1]]
  svdf_clump = data.frame(sv_list)
  svdf = data.frame(do.call('rbind', strsplit(as.character(svdf_clump$sv_list),'_',fixed=TRUE)))
  colnames(svdf) = c('start', 'end', 'sv')

  svdf$start = as.numeric(svdf$start)
  svdf$end = as.numeric(svdf$end)
  svdf$n = 1: dim(svdf)[1]

  svdfo = svdf
  for (i in 1:dim(svdf)[1]){
    if (svdf[i,'sv'] == 'del'){

      #start breakpoints after the del get a plus
      svdf[(svdf$start > svdf[i,'end']) & (svdf$n > i),'start'] =
        svdf[(svdf$start > svdf[i,'end']) & (svdf$n > i),'start'] + (svdf[i,'end'] - svdf[i,'start'])

      #end breakpoints after the del get a plus
      svdf[(svdf$end > svdf[i,'end']) & (svdf$n > i),'end'] =
        svdf[(svdf$end > svdf[i,'end']) & (svdf$n > i),'end'] + (svdf[i,'end'] - svdf[i,'start'])

    } else if (svdf[i,'sv'] == 'dup'){

      # start breakpoints after the dup get a minus
      svdf[(svdf$start > svdf[i,'end']) & (svdf$n > i),'start'] =
        svdf[(svdf$start > svdf[i,'end']) & (svdf$n > i),'start'] - (svdf[i,'end'] - svdf[i,'start'])

      # end breakpoints after thhe dup get a minus
      svdf[(svdf$end > svdf[i,'end']) & (svdf$n > i),'end'] =
        svdf[(svdf$end > svdf[i,'end']) & (svdf$n > i),'end'] - (svdf[i,'end'] - svdf[i,'start'])
    }
  }

  dfi = svdf[svdf$sv == 'inv',]
  any_overlap = F
  for (inv_n in 1:dim(dfi)[1]){
    inv_interval = data.frame(start = min(dfi[inv_n, c('start','end')]), end = max((dfi[inv_n, c('start','end')])))
    start_nested = any((svdf$start > inv_interval$start) & (svdf$start < inv_interval$end))
    end_nested =   any((svdf$end > inv_interval$start) & (svdf$end < inv_interval$end))

    if (start_nested | end_nested){
      any_overlap = T
    }
  }

  # print(svdfo)
  # print(svdf)
  # print(any_overlap)
  return(any_overlap)
}
