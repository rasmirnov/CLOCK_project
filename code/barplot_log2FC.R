library(dplyr)
library(tidyverse)
library(ggplot2)
setwd('/storage1/fs1/martyomov/Active/collaborations/romans/Greg_wu/')

# 1) load data
meta <- read.csv('CLOCK/output/output/scRNA/clock_meta.csv')

# 2. isolate cells from CSF & slice useless cols
path_figs <- 'CLOCK/output/figures/scRNA/barplot_log2FC/'
tissue_keep <- 'CSF'
# contrast_time <- 'W5/10'
contrast_time <- 'YR1'
meta_f <- meta %>% 
  filter(tissue %in% tissue_keep &
         time %in% c('BSE', contrast_time)) %>%  #'W05', 'W10', contrast_time
  # mutate(time = recode(time, 'W05' = 'W5_10', 'W10' = 'W5_10')) %>%
  select(time, cell_type, sample) 

# 3. calculate #cells per sample
num_tot_cells <- meta_f %>% group_by(sample) %>% 
  summarise(total_cells = n()) 

###----------------------------------------------------------------------------------------------
# 4.1. create list for each disease_group with unique samples: 
sampl_keep <- sapply(split(meta_f$sample, meta_f$time), unique)

# 4.2. count freq for all possible combination
df <- as.data.frame(table(meta_f$time,
                          meta_f$cell_type,
                          meta_f$sample))

# 4.3. rename colnames
df <- rename(df, time = Var1, cell_type = Var2,
             sample = Var3, cells_per_sample = Freq)

# 4.4. Convert disease_group to character vector
df$time <- as.character(df$time)

# 4.5.Filter rows based on sampl_keep
# exclude rows if sample was not initially presented in a disease_group
cell_type_counts <- subset(df, time %in% names(sampl_keep) & 
                             unlist(Map(function(x, y) x %in% y,
                                        df$sample, sampl_keep[df$time])))
###-----------------------------------------------------------------------------------------------

# 5. map sum column to adata_meta_per:
cell_type_counts$total_cells <- with(cell_type_counts,
                                     num_tot_cells$total_cells[match(sample, num_tot_cells$sample)])

# 6. calculate percent & rename 2 values
cell_type_counts <- cell_type_counts %>% 
  mutate(percent = cells_per_sample/total_cells * 100) 

# 7. calculate average percentage of each cell time in each group (e.g. BSE vs W5_10)
mean_cell_prop <- cell_type_counts %>% 
  group_by(time, cell_type) %>% 
  summarise(mean_prop = mean(percent))

# 8. compute log2(fold change)
fold_change_df <- mean_cell_prop %>%
  pivot_wider(names_from = 'time', values_from = 'mean_prop') %>% 
  mutate(Fold_Change = YR1/BSE) %>%   #W5_10
  mutate(log2FC = log2(Fold_Change)) %>% 
  subset(is.finite(log2FC))      # filter out Infinite values

# 9. filter out abs(log2FC) >= 0.5
fold_change_df <- fold_change_df %>% 
  filter(abs(log2FC) >= 0.5)

# 10. to avoid shift in colors when some variables are missing
fold_change_df$cell_type <- factor(fold_change_df$cell_type,
                                      levels=c('B Atypical', 'B Naive', 'B Switched', 'B Transitional',
                                               'CD14 Mono', 'CD16 Mono', 'CD4 Naive', 'CD4 Temra', 'CD4 Th1',
                                               'CD4 Th17', 'CD4 Th22', 'CD4 Treg memory', 'CD4 Treg naive',
                                               'CD4 exhausted', 'CD8 NKT-like', 'CD8 Naive', 'CD8 Tcm CCR4+',
                                               'CD8 Tcm CCR4-', 'CD8 Tem GZMB+', 'CD8 Tem GZMK+', 'CD8 Temra', 
                                               'CD8 Trm', 'ILC', 'MAIT', 
                                               'Microglia-like', 
                                               'NK CD56bright', 'NK CD56dim', 'Plasmablast', 
                                               'Proliferative T/NK', 'cDC', 'gdT2', 'pDC'))
my_palette <- c('#96705b', '#b69121', '#b76935', '#684756',   # B Transitional
                '#cdb4db', '#bcb6ff', '#d60000', '#f9cb9cff', '#e36414',
                '#ff9f1c', '#ffbf69','#edb8b8', '#ff7266', 
                '#fdf490',  '#74d3ae', '#a3b18a', '#3a5a40',
                '#588157',  '#007f5f', '#55a630', '#a5be00',
                '#bce784',  '#ff7ed1', '#00acc6', 
                '#d2d4c8',                               #MG
                '#eb0077', '#bf03b8', '#e7bc91',
                '#9651c4', '#595959', '#0774d8', '#a5a5a5')
names(my_palette) <- levels(fold_change_df$cell_type)

# 11. plot as a barplot
ggplot(fold_change_df, aes(x = reorder(cell_type, -log2FC), 
                           y = log2FC, fill = cell_type)) +
  geom_bar(stat = "identity", color='black') +
  theme_light() +
  theme(axis.title.x = element_blank(),
        legend.position="none",                     # alternative: guides(fill = "none")
        plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  scale_fill_manual(values = my_palette) +
  ggtitle(paste0(tissue_keep, ": BSE vs ", contrast_time)) +
  ylab(paste0("Log2 fold change (", contrast_time, " / BSE)")) 

# 12. save fig
ggsave(paste0(path_figs, tissue_keep, "_BSE_W5_10", ".pdf"),   #W5_10, contrast_time
       width = 6, height = 5)















