library(dplyr)
library(ggplot2)
library(readr)
library(tidyverse)

# 1. upload data
path <- '/storage1/fs1/martyomov/Active/collaborations/romans/Greg_wu'
atlas <- read.csv(paste0(path, "/CLOCK/data/single_cell/meta.csv"),
                       row.names = 1)
adata_meta <- atlas
tissue_v <- 'PBMC'

# 2. slice useless cols
adata_meta <- adata_meta %>%
  filter(tissue %in% tissue_v & time %in% c('BSE', 'W5', 'W10', 'YR1', 'YR2')) %>%  #, 'YR2'
  select(cell_type, time)

# 2.2. merge values from cols YR1 & YR2 to YR1_YR2
adata_meta$time <- sub("^(W).*", "W5_W10", adata_meta$time)
# adata_meta$time <- sub("^(Y).*", "YR1_YR2", adata_meta$time)

# 3. calculate percentages for each time point
adata <- adata_meta %>% group_by(time, cell_type) %>%    
  summarise(per_cluster = n()) %>%
  mutate(sum_cells = sum(per_cluster)) %>%
  mutate(percent = per_cluster/sum_cells*100)

# CHECK
adata %>% group_by(time) %>% 
  summarise(sum_perc = sum(percent))

# 4. fix the order of labels
adata$cell_type = factor(adata$cell_type,
                         levels=c('CD4 Naive', 'CD4 Treg naive', 'CD4 Treg memory', 'CD4 Th1',
                                  'CD4 Th17', 'CD4 Th22', 'CD4 Temra', 'CD4 exhausted',
                                  'CD8 Naive', 'CD8 Tcm CCR4-', 'CD8 Tcm CCR4+', 'CD8 Trm',
                                  'CD8 Tem GZMK+', 'CD8 Tem GZMB+', 'CD8 Temra', 'CD8 NKT-like',
                                  'gdT2', 'MAIT', 'Proliferative T/NK', 'NK CD56dim',
                                  'NK CD56bright', 'ILC', 'CD14 Mono', 'CD16 Mono',
                                  'cDC', 'pDC', 
                                  #'Microglia-like', 
                                  'B Naive','B Transitional', 'B Atypical', 'B Switched', 'Plasmablast'))

my_palette <- c('#d60000',  '#ff7266', '#edb8b8', '#e36414',
                '#ff9f1c', '#ffbf69', '#f9cb9cff', '#fdf490', 
                '#a3b18a', '#588157', '#3a5a40', '#bce784',
                '#55a630', '#007f5f', '#a5be00', '#74d3ae',
                '#0774d8', '#00acc6','#9651c4', '#bf03b8',  
                '#eb0077', '#ff7ed1',  '#cdb4db', '#bcb6ff', 
                '#595959', '#a5a5a5', 
                #'#d2d4c8', 
                '#b69121', '#684756', '#96705b', '#b76935', '#e7bc91')

# 5.1. plot stacked barplot
ggplot(adata, aes(x = time, y = percent, 
                  fill = as.factor(cell_type), label = cell_type))+    
  geom_bar(stat = "identity") +
  scale_fill_manual(name = 'Cell type', values = my_palette) +
  # scale_fill_discrete(name = 'Cell type') + 
  ggtitle(tissue_v) +
  ylab("Percentage of cells") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        axis.title.x=element_blank(),
        legend.text=element_text(size=12),
        legend.title=element_text(size=13))  
# facet_grid(rows = vars(time))

ggsave(paste0(path, "/CLOCK/output/figures/scRNA/freq_stacked_barpl/", tissue_v, "_all_sbarplot.png"), 
       dpi = 300, 
       width = 9, 
       height = 6)

ggsave(paste0(path, "/CLOCK/output/figures/scRNA/freq_stacked_barpl/", tissue_v, "_sbarplot.pdf"), 
       dpi = 300, 
       width = 9, 
       height = 6)

