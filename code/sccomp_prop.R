library(dplyr)
library(ggplot2)
library(forcats)
library(tidyr)
library(sccomp)
library(cmdstanr)
setwd('/storage1/fs1/martyomov/Active/collaborations/romans/Greg_wu/')


# 1) read my data
meta <- read.csv('CLOCK/output/output/scRNA/clock_meta.csv')

# 2. isolate cells from CSF & slice useless cols
path_figs <- 'CLOCK/output/figures/scRNA/abundance_volcano/all_labels/sccomp/'
tissue_keep <- 'PBMC'
# contrast_time <- 'W5/10'
contrast_time <- 'YR2'

meta_f <- meta %>% 
  filter(tissue %in% tissue_keep &
           time %in% c('BSE', contrast_time)) %>%  #'W05', 'W10', contrast_time
  # mutate(time = recode(time, 'W05' = 'W5_10', 'W10' = 'W5_10')) %>%
  select(time, cell_type, sample) 

# 3.1. count freq for all possible combination
df <- as.data.frame(table(meta_f$sample,
                          meta_f$cell_type))

# 3.2. rename colnames
df <- rename(df, sample = Var1, cell_type = Var2, cells_per_sample = Freq)

# 3.3. create time column
df$time <- vapply(strsplit(as.character(df$sample), '_'), 
       function(x) x[2], character(1))

# 3.4. make sure your count df has (1)#sample (2)time(type) (3)cell_type(cell_group) (4)count
my_counts <- df
my_counts <- my_counts %>% 
  rename(count = cells_per_sample) #%>% 
  # mutate(time = recode(time, 'W05' = 'W5_10', 'W10' = 'W5_10'))

# 4) run the scCOMP model
sccomp_res = 
  my_counts |>
  sccomp_estimate( 
    formula_composition = ~ time, 
    .sample = sample,
    .cell_group = cell_type,
    .count = count, 
    cores = 10, verbose = FALSE
  ) |> 
  sccomp_remove_outliers(cores = 10, verbose = FALSE) |> # Optional
  sccomp_test()

sccomp_res

# # 4.2) try to plot mine
sccomp_res |>
  sccomp_boxplot(factor = "time") #significance_threshold = 0.1
ggsave('CLOCK/output/figures/scRNA/abundance_volcano/all_labels/sccomp_boxplot/PDF/PBMC_YR2.pdf')

# 5) leave only relevant rows (without intercept)
sccomp_res_f <- sccomp_res %>% 
  filter(factor %in% "time") %>% 
  select(cell_type, c_FDR)

###----------------------------------------------------------------------------------------------

### compute Fold change
# 6) calculate #cells per sample
num_tot_cells <- meta_f %>% group_by(sample) %>% 
  summarise(total_cells = n()) 

# 7.1) create list for each disease_group with unique samples: 
sampl_keep <- sapply(split(meta_f$sample, meta_f$time), unique)

# 7.2) count freq for all possible combination
df <- as.data.frame(table(meta_f$time,
                          meta_f$cell_type,
                          meta_f$sample))

# 7.3) rename colnames
df <- rename(df, time = Var1, cell_type = Var2,
             sample = Var3, cells_per_sample = Freq)

# 7.4) Convert disease_group to character vector
df$time <- as.character(df$time)

# 7.5) Filter rows based on sampl_keep
# exclude rows if sample was not initially presented in a disease_group
cell_type_counts <- subset(df, time %in% names(sampl_keep) & 
                             unlist(Map(function(x, y) x %in% y,
                                        df$sample, sampl_keep[df$time])))

# 8) map sum column to adata_meta_per:
cell_type_counts$total_cells <- with(cell_type_counts,
                                     num_tot_cells$total_cells[match(sample, num_tot_cells$sample)])

# 9) calculate percent & rename 2 values
cell_type_counts <- cell_type_counts %>% 
  mutate(percent = cells_per_sample/total_cells * 100) 

# 10) calculate average percentage of each cell time in each group (e.g. BSE vs W5_10)
mean_cell_prop <- cell_type_counts %>% 
  group_by(time, cell_type) %>% 
  summarise(mean_prop = mean(percent))

# 11) compute log2(fold change)
fold_change_df <- mean_cell_prop %>%
  pivot_wider(names_from = 'time', values_from = 'mean_prop') %>% 
  mutate(fold_change = YR1/BSE) %>% 
  select(cell_type, fold_change)

# 12) join sccomp_res_f & fold_change_df by a common column
res <- merge(sccomp_res_f, fold_change_df, by=c("cell_type"))

# 13) reorder rows according to cell type
# for CSF: all 32
order <- c('B Atypical', 'B Naive', 'B Switched', 'B Transitional',
           'CD14 Mono', 'CD16 Mono', 'CD4 Naive', 'CD4 Temra', 'CD4 Th1',
           'CD4 Th17', 'CD4 Th22', 'CD4 Treg memory', 'CD4 Treg naive',
           'CD4 exhausted', 'CD8 NKT-like', 'CD8 Naive', 'CD8 Tcm CCR4+',
           'CD8 Tcm CCR4-', 'CD8 Tem GZMB+', 'CD8 Tem GZMK+', 'CD8 Temra',
           'CD8 Trm', 
           'ILC', 'MAIT', 
           'Microglia-like', 
           'NK CD56bright', 'NK CD56dim', 'Plasmablast', 'Proliferative T/NK', 'cDC', 'gdT2',
           'pDC')
res <- res %>% 
  slice(match(order, cell_type))
# delete cell type name colm
res <- res %>% select(c_FDR, fold_change)
# transpose
res_t <- t(res)
# res_t <- res_t[,-4]  # for CSF W5_10 only (remove transitional B)
# res_t <- res_t[,-22]  # for PBMC W5_10 & YR1 only (remove CD8 Trm)

# 14) save res
write.csv(res_t, 'CLOCK/output/output/scRNA/sccomp/CSF_YR1.csv', row.names=FALSE)

x <- read.csv('CLOCK/output/output/scRNA/sccomp/CSF_YR1.csv')
cell_types <- c('B Atypical', 'B Naive', 'B Switched', 'B Transitional',
                         'CD14 Mono', 'CD16 Mono', 'CD4 Naive', 'CD4 Temra', 'CD4 Th1',
                         'CD4 Th17', 'CD4 Th22', 'CD4 Treg memory', 'CD4 Treg naive',
                         'CD4 exhausted', 'CD8 NKT-like', 'CD8 Naive', 'CD8 Tcm CCR4+',
                         'CD8 Tcm CCR4-', 'CD8 Tem GZMB+', 'CD8 Tem GZMK+', 'CD8 Temra',
                         'CD8 Trm', 
                         'ILC', 'MAIT', 
                         'Microglia-like', 
                         'NK CD56bright', 'NK CD56dim', 'Plasmablast', 'Proliferative T/NK', 'cDC', 'gdT2',
                         'pDC')
colnames(x) <- cell_types
x_t <- t(x) %>% as.data.frame() %>% 
  mutate(log2FC = log2(V2))









