# Whoeps, 6th Apr 2022
# Explore output file.

library(dplyr)
library(stringr)

make_dumbbell_plot <- function(res1samp){


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
    geom_text(aes(x=105, y=n_line, label= paste0(simple_sample, ' ', mut_max)), size=3, hjust=0) +
    theme_bw() +
    scale_x_continuous(breaks=seq(0,100,10), limits=c(0,200)) +
    labs(x = 'Sequence concordance with hg38',
         title = paste0(res1samp[1,'invname']))


  return(p)
}


res_link = 'data/calls_combined.tsv'
res = read.table(res_link, sep='\t', header=T)

# # Hardcoded info
# samples_in = 33
# invs_in = 148
#
# barplot(table(res$mut_max))
#
# n_res = dim(res)[1]
# n_ref = dim(res[res$mut_max=='ref',])[1]
# n
# # N mutations
# res_ref = res[res$mut_max=='ref',]
# res_mut = res[res$mut_max!='ref',]
# hist(res_ref$res_max)
#
# cured = res[res$res_ref < 99 & res$res_max > 99,]
#
# hist(cured$res_ref)
# 3591
#
# res_nobreak = res[(res$exceeds_x == F) &  (res$exceeds_y == F),]
#
# res_ref %>% group_by(start)

starts = unique(res$start)
#starts = starts[1:5]
for (start in starts){
  res1 = res[res$start == start,]
  p = make_dumbbell_plot(res1)
  ggsave(filename=paste0('plots/', start, '.pdf'), width = 7.5, height=7.5, device='pdf', units='in')
}
