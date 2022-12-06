

# Times.
time_table = '~/PhD/projects/nahrcall/ntk-postprocess/stats.txt'
stats = read.table(time_table)
colnames(stats) = c('minutes', 'interval','len')

stats = stats[stats$len > 0,]


library(ggplot2)

ggplot(stats) + geom_point(aes(y=-minutes, x=len)) + theme_bw() +
  scale_x_log10() + labs(x='Inversion length', y='Processing time [minutes]',
                         title='Processing times per interval per sample.')

