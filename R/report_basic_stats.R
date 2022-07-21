# Whoeps, 6th Apr 2022
# Explore output file.

library(dplyr)
library(stringr)
library(ggplot2)

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
    geom_point(data=res1samp[res1samp$exceeds_y==T,], aes(x=res_ref, y=n_line), color='red') +
    geom_text(aes(x=105, y=n_line, label= paste0(simple_sample, ' ', mut_max)), size=3, hjust=0) +
    theme_bw() +
    scale_x_continuous(breaks=seq(0,100,10), limits=c(0,200)) +
    labs(x = 'Sequence concordance with hg38',
         title = paste0(res1samp[1,'invname']))

  p
  return(p)
}


#res_link = 'data/calls_combined.tsv'
#res_link = 'data/below50_7th_apr.tsv'

library(ggplot2)
res_link = 'data/above50_8th_apr.tsv'
res_link = '~/PhD/projects/nahrcall/all_samps_may12.tsv'
res_link = '~/Desktop/pics_fullrun_jul11/all.tsv'
res_link = '/Users/hoeps/PhD/projects/nahrcall/ntk-postprocess/res_15_jul/all.tsv'
res_link = '/Users/hoeps/PhD/projects/nahrcall/ntk-postprocess/res_15_jul/all_plusnachzug.tsv'
res_link = '/Users/hoeps/PhD/projects/nahrcall/analyses_paper/cao_regions/res/all.tsv'

res = read.table(res_link, sep='\t', header=T)

res = res[res$exceeds_y == F,]
res$size_mut = lengths(regmatches(res$mut_max, gregexpr("\\+", res$mut_max))) + 1

ape_samples = c('gorGor6_hg38', 'panPan3_hg38', 'panTro6_hg38', 'ponAbe3_hg38', 'rheMac10_hg38')

res = res[!(res$sample %in% ape_samples),]
#res$cp = paste0("cp ~/cluster16/scratch/hoeps/res_apr27/", res$seqname, "-", res$start, "-", res$end, '/diff/pdf/', res$sample, '*.png ~/Desktop/cure/', res$seqname, "-", res$start, "-", res$end, "-", res$sample,"-mut_",res$size_mut, ".png")
#res$cp_pdf = paste0("cp ~/cluster16/scratch/hoeps/res_apr27/", res$seqname, "-", res$start, "-", res$end, '/diff/pdf/', res$sample, '*.pdf ~/Desktop/cure/', res$seqname, "-", res$start, "-", res$end, "-", res$sample,"-mut_",res$size_mut, ".pdf")
#res$cp_grids = paste0("cp -r ~/cluster16/scratch/hoeps/res_apr27/", res$seqname, "-", res$start, "-", res$end, '/diff/pdf/grid ~/Desktop/cure/', res$seqname, "-", res$start, "-", res$end,  "/")


barplot(table(res$n_res_max))



res$inv_involved = lengths(regmatches(res$mut_max, gregexpr("inv", res$mut_max))) > 0
res[res$mut_max == 'ref',]$size_mut = 0
res$cured = res$res_max > 98


ggplot(res, aes(size_mut)) + geom_bar(aes(fill=inv_involved)) +
    scale_x_binned() +
    scale_x_continuous(breaks=c(0:4), limits=c(-1,4), labels=c('none','1mut','2mut','3mut','')) +
    labs(x='Number of mutations', y='Number of loci * samples',
         title='Number of mutations')

ggplot(res, aes(size_mut)) + geom_bar(aes(fill=cured)) +
  scale_x_binned() +
  scale_x_continuous(breaks=c(0:4), limits=c(-1,4), labels=c('none','1mut','2mut','3mut','')) +
  labs(x='Number of mutations', y='Number of loci * samples',
       title='>99 sequence agreement post mutation')

rc = res[res$cured == T,]

res = res[order(res$seqname, res$start),]

rc = rc[order(rc$seqname, rc$start),]

res_cure_4 = res[((res$cured==T)&(res$size_mut==4)),]
res_cure_3 = res[((res$cured==T)&(res$size_mut==3)),]
res_cure_2 = res[((res$cured==T)&(res$size_mut==2)),]

table(res$cured, res$size_mut)


#numbers = data.frame(mut = res$mut_max, n = n_stuff)





n_res = dim(res)[1]
n_ref = dim(res[res$mut_max=='ref',])[1]

# N mutations
res_ref = res[res$mut_max=='ref',]
res_mut = res[res$mut_max!='ref',]
hist(res_ref$res_max)

cured = res[res$res_ref < 98 & res$res_max > 98,]

hist(cured$res_ref)


res_nobreak = res[(res$exceeds_x == F) &  (res$exceeds_y == F),]

res_ref %>% group_by(start)

intervals = unique(res[,c('seqname', 'start','end')])
n = 1
chr = intervals[n, 'seqname']
start = intervals[n, 'start']
end = intervals[n, 'end']


res$inv_involved = lengths(regmatches(res$mut_max, gregexpr("inv", res$mut_max))) > 0
res_cure = res[(res$cured==T) & (res$size_mut > 1),]
res_cure$onlydels = lengths(regmatches(res_cure$mut_max, gregexpr("del", res_cure$mut_max))) == res_cure$size_mut
res_cure$onlydups = lengths(regmatches(res_cure$mut_max, gregexpr("dup", res_cure$mut_max))) == res_cure$size_mut
res_cure_interest = res_cure[(res_cure$onlydels + res_cure$onlydups == 0) & (res_cure$inv_involved == T),]
# intervals = unique(res[,c('seqname', 'start','end')])
# n_intervals = dim(intervals)[1]
# #n_intervals = 1
# for (n in 1:n_intervals){
#
#   chr = intervals[n, 'seqname']
#   start = intervals[n, 'start']
#   end = intervals[n, 'end']
#   print('hi')
#   res1 = res[res$start == start,]
#   p = make_dumbbell_plot(res1)
#
#   ggsave(filename=paste0('plots_above50/',chr,'-', start,'-', end, '.pdf'), width = 7.5, height=7.5, device='pdf', units='in')
# }

#write.table(file='~/Desktop/copy2.sh', res_cure_2$cp, quote=F, row.names=F, col.names=F)
#write.table(file='~/Desktop/copy3.sh', res_cure_3$cp, quote=F, row.names=F, col.names=F)
#write.table(file='~/Desktop/copy4.sh', res_cure_4$cp, quote=F, row.names=F, col.names=F)

#write.table(file='~/Desktop/copy2p.sh', res_cure_2$cp_pdf, quote=F, row.names=F, col.names=F)
#write.table(file='~/Desktop/copy3p.sh', res_cure_3$cp_pdf, quote=F, row.names=F, col.names=F)
#write.table(file='~/Desktop/copy4p.sh', res_cure_4$cp_pdf, quote=F, row.names=F, col.names=F)

#write.table(file='~/Desktop/copy2g.sh', res_cure_2$cp_grids, quote=F, row.names=F, col.names=F)
#write.table(file='~/Desktop/copy3g.sh', res_cure_3$cp_grids, quote=F, row.names=F, col.names=F)
#write.table(file='~/Desktop/copy4g.sh', res_cure_4$cp_grids, quote=F, row.names=F, col.names=F)

#
# save2 = res_cure_2 %>% group_by(start) %>% slice(1)
# save3 = res_cure_3 %>% group_by(start) %>% slice(1)
# save4 = res_cure_4 %>% group_by(start) %>% slice(1)
#
# write.table(rc %>% group_by(start) %>% slice(1), file='~/Desktop/cure/all.tsv', sep='\t', quote=F, row.names=F, col.names=F)
# write.table(save2[order(save2$seqname, save2$start),], file='~/Desktop/cure/rc2', sep='\t', quote=F, row.names=F, col.names=T)
# write.table(save3[order(save3$seqname, save3$start),], file='~/Desktop/cure/rc3', sep='\t', quote=F, row.names=F, col.names=T)
# write.table(save4[order(save4$seqname, save4$start),], file='~/Desktop/cure/rc4', sep='\t', quote=F, row.names=F, col.names=T)

uniqs = res_cure_interest %>% group_by(start) %>% slice(1)
uniqs = uniqs[uniqs$end - uniqs$start > 20000,]
a = uniqs[c('seqname','start','end')]
write.table(a,'~/Desktop/pics_fullrun_jul11/cp.sh', sep='-', quote=F, row.names=F, col.names=F)
