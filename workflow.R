setwd("/home/lty/PDE_Func_Est")

source('simu_begin.R')
source('simu_rand.R')

rm(list=ls())
gc()

source('real_prep.R')
save.image('real_prep.RData')
load('real_prep.RData')

source('real_cut.R')
save(list = c('ini.coef','ZTZ','ZTu','k_list','BIC_cut','idx_K','ucut_frap'), file = 'real_cut.RData')
load('real_cut.RData')

source('real_pen.R')
save(list = c('b.coef','bTb','bTu','pen_bd','pen_pde','pen','tune','GCV_pen','idx_tune','upen_frap'), file = 'real_pen.RData')
load('real_pen.RData')

source('real_sum.R')

