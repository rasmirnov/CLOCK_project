library(dplyr)
library(ggplot2)
library(tidyverse)

# 1.1 upload data
path <- '/storage1/fs1/martyomov/Active/collaborations/romans/Greg_wu/'
atlas <- read.csv(paste0(path, "/CLOCK/data/meta.csv"),
                  row.names = 1)
adata_meta <- atlas
organ <- 'CSF'

# 2. slice useless cols & filter out samples with low #cells
adata_meta <- adata_meta %>%
  # filter((tissue %in% organ) & (time %in% c("BSE", "YR2"))) %>% 
  filter(tissue %in% organ) %>% 
  select(cell_type, sample,  time)

# 2.2. check samples with low #cells
tot_per_sample <- adata_meta %>% 
  group_by(sample) %>% 
  summarise(total_cells = n())  
tot_per_sample

# 2.3 slice useless cols & filter out samples with low #cells
adata_meta <- adata_meta %>%
  filter(!sample %in% c('009_BSE_CSF', '011_BSE_CSF'))

# 2.4. merge YR1 & YR2 as well as W10 & W5
adata_meta <- adata_meta %>% 
  mutate(time = recode(time, 'YR1' = 'YR1_YR2', 'YR2' = 'YR1_YR2',
                       'W05' = 'W5_W10', 'W10' = 'W5_W10'))

# 3. count total #cells per sample
tot_per_time <- adata_meta %>% 
  group_by(time) %>% 
  summarise(total_cells = n())    

# 4. count #cells per cell_type and Time group
adata <- adata_meta %>% 
  mutate(cell_type=as.factor(cell_type), time=as.factor(time)) %>% 
  count(cell_type, time, .drop = F, name = 'per_cell_type_time')

# 5. map total_cells column to adata dataset:
# adata$tot_per_time <- with(adata, tot_per_time$total_cells[match(time, tot_per_time$time)])
adata <- merge(adata, tot_per_time, by=c("time"))

# 6. count # unique samples per "time" group
n_samples <- adata_meta %>% 
  group_by(time) %>% 
  summarise(n_samples_per_time = n_distinct(sample))  

# 7. map n_samples_per_time column to adata dataset:
# adata$n_samples_per_time <- with(adata, n_samples$n_samples_per_time[match(time, n_samples$time)])
adata <- merge(adata, n_samples, by=c("time"))

# 8. calculate "fraction" of cells per cell type X across "time" group
adata <- adata %>% 
  group_by(cell_type, time) %>% 
  mutate(fraction = per_cell_type_time/(total_cells * n_samples_per_time))

# 9. calculate "percent" of cells per cell type X across "time" group
adata <- adata %>% 
  group_by(cell_type) %>% 
  mutate(sum_frac = sum(fraction)) %>%
  mutate(percent = fraction/sum_frac*100)

# 10. fix the order of labels
adata$cell_type = factor(adata$cell_type,
                         levels=c('CD4 Naive', 'CD4 Treg naive', 'CD4 Treg memory', 'CD4 Th1', 'CD4 Th17', 'CD4 Th22', 'CD4 Temra',
                                  'CD4 exhausted', 'CD8 Naive', 'CD8 Tcm CCR4-', 'CD8 Tem GZMK+', 'CD8 Tem GZMB+', 'gdT1', 'gdT2', 
                                  'MAIT', 'Proliferative T/NK', 'NK CD56dim', 'NK CD56bright', 'ILC', 'CD14 Mono', 'CD16 Mono', 'cDC', 'pDC', 
                                  'Microglia-like', 'B Naive', 'B Transitional', 'B Atypical', 'B Switched', 'Plasmablast'))

# FOR CSF
# 'CD4 Naive', 'CD4 Treg naive', 'CD4 Treg memory', 'CD4 Th1', 'CD4 Th17', 'CD4 Th22', 'CD4 Temra',
# 'CD4 exhausted', 'CD8 Naive', 'CD8 Tcm CCR4-', 'CD8 Tem GZMK+', 'CD8 Tem GZMB+', 'gdT1', 'gdT2', 
# 'MAIT', 'Proliferative T/NK', 'NK CD56dim', 'NK CD56bright', 'ILC', 'CD14 Mono', 'CD16 Mono', 'cDC', 'pDC', 
# 'Microglia-like', 'B Naive', 'B Transitional', 'B Atypical', 'B Switched', 'Plasmablast'

# FOR PBMC
# my_palette <- c('#d60000',  '#ff7266', '#edb8b8', '#e36414', '#ff9f1c', '#ffbf69', '#f9cb9cff', '#fdf490', 
#                 '#018700', '#004b00',  '#5d7e66', '#95b577', '#0000dd', 
#                 '#0774d8', '#00acc6','#9651c4', '#bf03b8',  
#                 '#eb0077', '#ff7ed1',  '#cdb4db', '#bcb6ff', '#595959', '#a5a5a5', 
#                 '#b69121', '#684756', '#96705b', '#b76935', '#e7bc91')
# 11. define my_palette and fix the order of colors
my_palette <- c('#d60000',  '#ff7266', '#edb8b8', '#e36414', '#ff9f1c', '#ffbf69', '#f9cb9cff', '#fdf490', 
                '#018700', '#004b00',  '#5d7e66', '#95b577', '#0000dd', 
                '#0774d8', '#00acc6','#9651c4', '#bf03b8',  
                '#eb0077', '#ff7ed1',  '#cdb4db', '#bcb6ff', '#595959', '#a5a5a5', 
                '#d2d4c8', '#b69121', '#684756', '#96705b', '#b76935', '#e7bc91')

# FOR CSF

# FOR PBMC

# 12. plot dotplot
ggplot(adata, aes(x = time, y = cell_type, 
                  size = percent, colour = cell_type))+  
  geom_point() +
  # scale_color_discrete(name = "Cell type") +
  scale_color_manual(values = my_palette,
                     name = 'Cell type') +
  ggtitle(organ) +
  guides(colour = "none") +    # remove a color-legend
  xlab("Time of treatment") +
  ylab("Percentage of cells") +
  theme_bw() +
  theme(panel.grid = element_line(color = rgb(235, 235, 235, 150, maxColorValue = 500)),  #col_grid
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle=45, hjust=1),
        axis.text=element_text(size=10),
        axis.title=element_text(size=12),
        legend.text=element_text(size=10),
        legend.title=element_text(size=12)) +
  scale_y_discrete(limits=rev)

ggsave(paste0("/storage1/fs1/martyomov/Active/collaborations/romans/Greg_wu/CLOCK/output/figures/",
              'png/', organ, "_dot_plot_all.png"), 
       dpi = 300, 
       width = 4, 
       height = 7)

ggsave(paste0("~/Downloads/MD_proj/CLOCK_proj/figures/", organ, "_dot_plot.pdf"), 
       dpi = 300, 
       width = 4, 
       height = 6)