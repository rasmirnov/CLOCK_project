library(ggplot2)
library(reticulate)
use_python('/storage1/fs1/martyomov/Active/collaborations/romans/Greg_wu/CSF_atlas/code/scanpy/figs_de/r-reticulate/bin/python3.8')
library(Seurat)
library(anndata)
library(dplyr)
library(tidyverse)
library(viridis)
library(data.table)

# 1. set env for loading the data
path <- "/storage1/fs1/martyomov/Active/collaborations/romans/Greg_wu/CLOCK/"
setwd(paste0(path, 'output/output/'))

# 2. upload raw h5ad obj
adata <- read_h5ad("combined_annot.h5ad")

# 3. create seurat obj with raw counts
seurat_data <- CreateSeuratObject(counts = t(as.matrix(adata$layers[['counts']])),      
                                  meta.data = adata$obs)
# 4. check raw counts
seurat_data@assays$RNA@counts %>% max()

# 5. merge W5 & W10 or YR1 & YR2
seurat_data@meta.data <- seurat_data@meta.data %>% 
  mutate(cell_type = recode(cell_type, 'Proliferative T/NK' = 'Proliferative T&NK'))

# 6.1. subset based on matrix design
tissue_s <- 'CSF'
seurat_subset <- subset(seurat_data, tissue %in% tissue_s)
seurat_subset <- subset(seurat_subset, time %in% c('BSE', 'W05', 'W10', 'YR1', 'YR2'))
cell_type_list <- c('CD8 Naive', 'CD8 Tcm CCR4-', 'CD8 Tem GZMK+', 'CD8 Tem GZMB+', 'gdT1', 'gdT2', 
                    'MAIT') # remove only for CSF-W5_10: 'B Transitional', 
seurat_subset <- subset(seurat_subset, cell_type %in% cell_type_list)

# 6.2. add extra column with cell_type_time
# seurat_subset@meta.data$cell_type.time <- paste0(seurat_subset@meta.data$cell_type, "_", seurat_subset@meta.data$time)

# 7. make averaged sum of aggregated raw counts
pbulk <- AggregateExpression(
  object = seurat_subset, assays = 'RNA',
  features = NULL,  return.seurat = TRUE,
  group.by = c('cell_type', 'sample'),
  add.ident = NULL, slot = "counts", verbose = TRUE)

# 4. normilize (for DE)
pbulk <- NormalizeData(pbulk, normalization.method = "LogNormalize",
                       scale.factor = 10000)

# 5. scale (for heatmap)
all.genes <- rownames(pbulk)
pbulk <- ScaleData(pbulk, features = all.genes)

# 8.1 check
# pbulk@meta.data <- seurat_subset@meta.data
# levels(pbulk)
# Idents(pbulk)
# pbulk@assays$RNA@counts %>% dim()
# pbulk@assays$RNA@counts %>% max()

# 8.2. modify data
pbulk@meta.data$cell_type.subject.time.tissue <- rownames(pbulk@meta.data)
pbulk@meta.data <- pbulk@meta.data %>%
  mutate(time = sapply(strsplit(as.character(cell_type.subject.time.tissue), split = "_"), `[`, 3)) %>%
  mutate(time = recode(time, 'W05' = 'W5_10', 'W10' = 'W5_10', 'YR1' = 'YR1_2', 'YR2' = 'YR1_2')) %>% 
  rename(cell_type = orig.ident) %>% 
  mutate(cell_type.time = paste0(cell_type, '_', time))

pbulk@meta.data$cell_type.time %>% unique()
###-------------------------------------------Using pheatmap----------------------------------------------------
## using pheatmap for multiply bar legend (after pbulk & DEA steps)
library(ComplexHeatmap)
library(pheatmap)
library(RColorBrewer)

### 9. load all DEGs per tissue-treatment group 
degs_w <- read.csv(paste0('DEA/', tissue_s, '_W5_10/combined.csv'))
degs_yr <- read.csv(paste0('DEA/', tissue_s, '_YR1_2/combined.csv'))
all_degs <- bind_rows(degs_w, degs_yr)

ms_risk_genes <- read.csv('MS_risk_genes/GWAS_MS_risk_genes.csv')

# 10. refine metadata
all_degs <- all_degs %>% 
  rename(gene = V1) %>%
  subset(select = -X)

# 11. leave only DE MS risk genes for a specific cluster
sub_degs <- all_degs %>% 
  filter(gene %in% ms_risk_genes$risk_genes) %>% 
  filter(p_val < 0.05) %>%      
  filter(avg_log2FC > 0.5 | avg_log2FC < -0.5) %>% 
  filter(cell_type %in% cell_type_list) %>% 
  distinct(gene, .keep_all = TRUE)   # remove duplicated genes

# 12. order labels
# levels(pbulk) <- c("B Naive_BSE",       "B Naive_W5_10", "B Naive_YR1_2",
#                    "B Atypical_BSE", "B Atypical_W5_10", "B Atypical_YR1_2",
#                    "B Switched_BSE", "B Switched_W5_10", "B Switched_YR1_2",
#                    "Plasmablast_BSE", "Plasmablast_W5_10", "Plasmablast_YR1_2")
pbulk$cell_type.time <- factor(pbulk$cell_type.time,
                               levels = c("CD8 Naive_BSE",      "CD8 Naive_W5_10",     "CD8 Naive_YR1_2",    
                                 "CD8 Tcm CCR4-_BSE",  "CD8 Tcm CCR4-_W5_10", "CD8 Tcm CCR4-_YR1_2",
                                 "CD8 Tem GZMB+_BSE",  "CD8 Tem GZMB+_W5_10", "CD8 Tem GZMB+_YR1_2",
                                 "CD8 Tem GZMK+_BSE", "CD8 Tem GZMK+_W5_10", "CD8 Tem GZMK+_YR1_2",
                                 "gdT1_BSE",           "gdT1_W5_10",          "gdT1_YR1_2",         
                                 "gdT2_BSE",           "gdT2_W5_10",          "gdT2_YR1_2",         
                                 "MAIT_BSE",           "MAIT_W5_10",          "MAIT_YR1_2"))
levels(pbulk$cell_type.time)


# 13. extract gene expression matrix selecting rows that correspond to marker genes
DefaultAssay(pbulk) = "RNA"
datamat = pbulk$RNA@scale.data[which(rownames(pbulk) %in% sub_degs$gene), ]    
datamat = datamat[match(sub_degs$gene, rownames(datamat)), ]
datamat = datamat[, order(pbulk$cell_type.time)]
celltype = pbulk$cell_type.time
celltype = celltype[order(pbulk$cell_type.time)]

annotation_col = data.frame(cell_type_time = factor(celltype))      
rownames(annotation_col) = colnames(datamat)

# 14. create 2 separate colms: cell_type & organ 
annotation_col <- annotation_col %>% 
  mutate('Cell type' = sapply(cell_type_time, function(x) {gsub("_.*", "", x)})) %>% 
  mutate(Treatment = ifelse(grepl('BSE', cell_type_time), 'BSE',
                            ifelse(grepl('W5_10', cell_type_time), 'W5_10', 'YR1_2'))) %>%
  select('Cell type', Treatment)

# 15. reorder columns to make "cell_type" on the top
annotation_col <- annotation_col[ ,c(2, 1)]

# 16. clip heatmap values to plus or minus 2.5, to ensure that color contrast is high
datamat[datamat <= -2.3] = -2.3
datamat[datamat >= 2.3] = 2.3

# 17. choose colors
mycolors <- c('#018700', '#004b00',  '#5d7e66', '#95b577', '#0000dd', '#0774d8',
              '#00acc6')
names(mycolors) <- unique(annotation_col$'Cell type')
ann_colors <- list('Cell type' = mycolors,
                   Treatment = c(BSE = '#cc92c2', W5_10 = '#aa7bc3', YR1_2 = '#944bbb'))   # for YR1_2: #944bbb
# 18. order
annotation_col$'Cell type'  = factor(annotation_col$'Cell type',
                                     levels = cell_type_list)

# 19. plot
mapal <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256)
p <- ComplexHeatmap::pheatmap(datamat,         
              color = rev(mapal), 
              border_color = "grey60",
              name = "Expression",         # scale title
              annotation_col = annotation_col,
              annotation_colors = ann_colors,
              show_colnames = F, cluster_rows = F, cluster_cols = F,
              annotation_legend = T, 
              column_split = annotation_col$'Cell type', #gaps_col = T,
              column_title_gp = gpar(fontsize = 10), column_title_rot = 45,
              fontsize = 10, fontsize_row = 5)  

# 20. save
path_2 <- '/storage1/fs1/martyomov/Active/collaborations/romans/Greg_wu/CLOCK/output/figures/pdf/DEGs_risk_genes/pbulk/'
# install.packages('tidyHeatmap')
# library(tidyHeatmap)
save_pdf(p, paste0(path_2, "heatmap_", tissue_s, "_CD8.pdf"),
  width = 9, height = 7)

# # 19. define the function
# save_pheatmap_pdf <- function(x, filename, width=7, height=6) {
#   stopifnot(!missing(x))
#   stopifnot(!missing(filename))
#   pdf(filename, width=width, height=height)
#   grid::grid.newpage()
#   grid::grid.draw(x$gtable)
#   dev.off() }
###--------------------------------------------All labels dictionary------------------------------
# "CD4 Naive_BSE",         "CD4 Naive_W5_10",       "CD4 Naive_YR1_2",
# "CD4 Treg naive_BSE",    "CD4 Treg naive_W5_10",  "CD4 Treg naive_YR1_2",
# "CD4 Treg memory_BSE",   "CD4 Treg memory_W5_10", "CD4 Treg memory_YR1_2",
# "CD4 Th1_BSE",           "CD4 Th1_W5_10",         "CD4 Th1_YR1_2",        
# "CD4 Th17_BSE",          "CD4 Th17_W5_10",        "CD4 Th17_YR1_2",       
# "CD4 Th22_BSE",          "CD4 Th22_W5_10",        "CD4 Th22_YR1_2",       
# "CD4 Temra_BSE",         "CD4 Temra_W5_10",       "CD4 Temra_YR1_2",
# "CD4 exhausted_BSE",     "CD4 exhausted_W5_10",   "CD4 exhausted_YR1_2"

# "CD8 Naive_BSE",      "CD8 Naive_W5_10",     "CD8 Naive_YR1_2",    
# "CD8 Tcm CCR4-_BSE",  "CD8 Tcm CCR4-_W5_10", "CD8 Tcm CCR4-_YR1_2",
# "CD8 Tem GZMB+_BSE",  "CD8 Tem GZMB+_W5_10", "CD8 Tem GZMB+_YR1_2",
# "CD8 Tem GZMK+_BSE", "CD8 Tem GZMK+_W5_10", "CD8 Tem GZMK+_YR1_2",
# "gdT1_BSE",           "gdT1_W5_10",          "gdT1_YR1_2",         
# "gdT2_BSE",           "gdT2_W5_10",          "gdT2_YR1_2",         
# "MAIT_BSE",           "MAIT_W5_10",          "MAIT_YR1_2"

# "Proliferative T&NK_BSE",       "Proliferative T&NK_W5_10", "Proliferative T&NK_YR1_2",
# "NK CD56dim_BSE", "NK CD56dim_W5_10", "NK CD56dim_YR1_2",
# "NK CD56bright_BSE", "NK CD56bright_W5_10", "NK CD56bright_YR1_2",
# "ILC_BSE", "ILC_W5_10", "ILC_YR1_2"

# "CD14 Mono_BSE",       "CD14 Mono_W5_10", "CD14 Mono_YR1_2",
# "CD16 Mono_BSE", "CD16 Mono_W5_10", "CD16 Mono_YR1_2",
# "cDC_BSE", "cDC_W5_10", "cDC_YR1_2",
# "pDC_BSE", "pDC_W5_10", "pDC_YR1_2",
# "Microglia-like_BSE", "Microglia-like_W5_10", "Microglia-like_YR1_2"

# "B Naive_BSE",       "B Naive_W5_10", "B Naive_YR1_2",
# "B Atypical_BSE", "B Atypical_W5_10", "B Atypical_YR1_2",
# "B Switched_BSE", "B Switched_W5_10", "B Switched_YR1_2",
# "Plasmablast_BSE", "Plasmablast_W5_10", "Plasmablast_YR1_2"