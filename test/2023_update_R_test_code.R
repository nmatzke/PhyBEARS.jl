#######################################################
# Update to R code after 2023 diversitree (apparently)
#######################################################

lq_including_root =  = attr(res,"intermediates")$lq
tmp_lq = lq_including_root
tmp_lq[rootnode] = 0
lq = t(tmp_lq)
lq
