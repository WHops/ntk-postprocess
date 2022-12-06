

add_n_ape_resolved <- function(res_f, res_th = 98){

res_ape = res_f[res_f$sample %in% ape_samples,]
res_ape$ape_resolve = F
res_ape[res_ape$res_max > res_th,]$ape_resolve = T
n_resolve = res_ape %>% group_by(start) %>% mutate(sum(ape_resolve==T)) %>% slice(1)
n_resolve_sliced = res_ape %>% group_by(start) %>% mutate(n_ape_resolved = sum(ape_resolve==T)) %>% slice(1)

n_resolve_sliced_condensed = n_resolve_sliced[,c('start','n_ape_resolved')]
res_f_apes = left_join(res_f, n_resolve_sliced_condensed, by='start')

return(res_f_apes)
}
