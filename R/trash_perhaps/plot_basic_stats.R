# Plots
library(reshape2)
library(ggplot2)
library(pheatmap)

ancestry_file = '/Users/hoeps/Desktop/pics_20_selected/ancestries.sh'
anc = read.table(ancestry_file, sep='\t', col.names = c('simplesample','ANC'))



# Experiment
plot_one_locus_v2(res, anc, start_coord, solve_th)

for (start_coord in unique(res$start)){
  print(start_coord)
#start_coord = '25682654'
solve_th = 98
plot1(res, start_coord, solve_th)
}

