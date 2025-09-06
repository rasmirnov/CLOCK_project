library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyverse)

# 1. upload data
atlas_meta <- read.csv("~/Downloads/MD_proj/CLOCK_proj/data/clock_meta.csv",
                       row.names = 1)
adata_meta <- atlas_meta
organ_v <- 'CSF'

# add 2 more columns for simplicity
adata_meta$time <- adata_meta$Drug
adata_meta$organ <- sapply(adata_meta$Drug_Tissue,
                            function(x) {str_extract(x, "[^_]*")})

# 2. slice useless cols
adata_meta <- adata_meta %>%
  filter(organ %in% organ_v & time %in% c('BSE', 'YR1')) %>% 
  select(cell_type, time)

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
                         levels=c('CD4 Naive', 'CD4 Treg naive', 'CD4 Treg memory', 'CD4 Th1', 'CD4 Th17', 'CD4 Th22', 'CD4 Proliferating',
                                  'CD8 Naive', 'CD8 Tcm CCR4-', 'CD8 Tem GZMK+', 'CD8 Tem GZMB+', 'gdT1', 'gdT2', 'MAIT', 'NK CD56dim',
                                  'NK CD56bright', 'CD14 Mono', 'CD16 Mono', 'cDC1', 'cDC2', 'pDC', 'Microglia-like', 'B Naive',
                                  'B Switched', 'Plasmablast', 'MAST', 'RBC', 'Platelet'))

my_palette <- c('#d60000',  '#ff7266', '#edb8b8', '#ffa52f', '#f1c232ff', '#f9cb9cff', '#fdf490',  
                '#018700', '#004b00',  '#5d7e66', '#95b577', '#0000dd', 
                '#0774d8', '#00acc6', '#30cdcd', '#9ae4ff', '#9651c4', '#bf03b8', 
                '#eb0077', '#ff7ed1',  '#cdb4db', '#bcb6ff', '#684756', '#96705b', '#e7bc91', 
                '#595959', '#a5a5a5', '#d2d4c8')
  
# 5.1. plot stacked barplot
ggplot(adata, aes(x = time, y = percent, 
                  fill = as.factor(cell_type), label = cell_type))+    
  geom_bar(stat = "identity") +
  scale_fill_manual(name = 'Cell type', values = my_palette) +
  # scale_fill_discrete(name = 'Cell type') + 
  ggtitle(organ_v) +
  xlab("Time of treatment") +
  ylab("Percentage of cells") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        legend.text=element_text(size=12),
        legend.title=element_text(size=13))  
  # facet_grid(rows = vars(time))


ggsave(paste0("~/Downloads/MD_proj/CLOCK_proj/figures/", organ_v, "_sbarplot.png"), 
       dpi = 300, 
       width = 9, 
       height = 6)

ggsave(paste0("~/Downloads/MD_proj/CLOCK_proj/figures/", organ_v, "_sbarplot.pdf"), 
       dpi = 300, 
       width = 9, 
       height = 6)

