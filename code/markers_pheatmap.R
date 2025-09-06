library(ggplot2)
library(Seurat)
library(anndata)
library(dplyr)
library(tidyverse)
library(viridis)
library(data.table)
library(ComplexHeatmap)
library(pheatmap)
library(RColorBrewer)
library(tidyHeatmap)

# 1) set env for loading the data & load the obj
path <- "/storage1/fs1/martyomov/Active/collaborations/romans/Greg_wu/CLOCK/output/output/"
setwd(path)
expData <- readRDS('scRNA/combined_annot.rds')

# 1) rename Proliferative T/NK cell type
expData_s <- expData
expData_s@meta.data <- expData_s@meta.data %>% 
  mutate(cell_type = recode(cell_type, 'Proliferative T/NK' = 'Proliferative T&NK')) 

# 1) select specific cell types
cell_type_list <- c('CD4 Naive', 'CD4 Treg naive', 'CD4 Treg memory', 'CD4 Th1', 
  'CD4 Th17', 'CD4 Th22', 'CD4 Temra', 'CD4 exhausted', 'CD8 Naive', 'CD8 Tcm CCR4-',
  'CD8 Tcm CCR4+', 'CD8 Trm', 'CD8 Tem GZMK+', 'CD8 Tem GZMB+', 'CD8 Temra', 'CD8 NKT-like',
  'gdT2', 'MAIT')
# 'Proliferative T&NK', 'NK CD56dim', 'NK CD56bright',
# 'ILC', 'CD14 Mono', 'CD16 Mono', 'cDC', 'pDC', 'Microglia-like', 
# 'B Transitional', 'B Naive', 'B Atypical', 'B Switched', 'Plasmablast'

# 'CD4 Naive', 'CD4 Treg naive', 'CD4 Treg memory', 'CD4 Th1', 
# 'CD4 Th17', 'CD4 Th22', 'CD4 Temra', 'CD4 exhausted', 'CD8 Naive', 'CD8 Tcm CCR4-',
# 'CD8 Tcm CCR4+', 'CD8 Trm', 'CD8 Tem GZMK+', 'CD8 Tem GZMB+', 'CD8 Temra', 'CD8 NKT-like',
# 'gdT2', 'MAIT'

# 1) subset based on cell types
expData_s <- subset(expData_s, cell_type %in% cell_type_list)

# 1) make averaged sum of aggregated raw counts
expData_s <- AverageExpression(
  object = expData_s,
  assays = 'RNA', features = NULL,
  return.seurat = TRUE, group.by = c('cell_type'),  
  add.ident = NULL, slot = "counts", verbose = TRUE)

# 1) check
expData_s@assays$RNA@counts %>% dim()
expData_s@assays$RNA@counts %>% max()

# 1) normilize (for heatmap)
expData_s <- NormalizeData(expData_s, normalization.method = "LogNormalize",
                       scale.factor = 10000)   #pbulk

# 1) modify data
expData_s@meta.data$cell_type <- rownames(expData_s@meta.data)

# 1) order labels
levels(expData_s) <- cell_type_list
expData_s$cell_type <- factor(expData_s$cell_type, levels = cell_type_list)
levels(expData_s$cell_type)

# check activation/inflamation markers
# FeaturePlot(expData, c('IL17A', 'IL17B', 'IL17C', 'IL17D', 'IL17E', 'IL17F', 'IFNG', 'TNF'))

# 1) extract gene expression matrix selecting rows that correspond to marker genes
markers_lst <- data.frame(gene = c('CD3E', 'CD4', 'CCR7', 
               'FOXP3', 'TCF7','HLA-DRB1', 'HLA-DRB5', 'FANK1',
               'KLRB1', 'GZMA', 'CXCR3','RORC', 'CCR10', 'CCR6',    #'SEMA3G',
               'GZMH', 'PRF1', 'GNLY', 'MYB', 'CTLA4', 'PDCD1',
               'CD8A', 'CD8B', 'CCR4', 'ITGA1', 'CXCR6', 'GZMK', 'EOMES', 'KLRK1',
               'DUSP2', 'GZMB', 'KLRC2', 'IKZF2', 'TYROBP',   # 'NCAM1',
               'TRDV2', 'TRAV1-2', 'SLC4A10', 'CEBPD', 'IFNG', 'IL17D'))
# 'CD3E', 'CD4', 'CCR7', 'FOXP3', 'TCF7','HLA-DRB1', 'HLA-DRB5', 'FANK1', 'SEMA3G','KLRB1', 'GZMA', 'CXCR3','RORC', 'CCR10', 'CCR6', 'GZMH', 'PRF1', 'GNLY', 'MYB', 'CTLA4', 'PDCD1',
# 'CD8A', 'CD8B', 'CCR4', 'ITGA1', 'CXCR6', 'GZMK', 'EOMES', 'KLRK1', 'DUSP2', 'GZMB', 'KLRC2', 'IKZF2', 'NCAM1', 'TYROBP', 'TRDV2', 'TRAV1-2', 'SLC4A10', 'CEBPD'

# 'CD3E','MKI67', 'TYMS','NCR1', 'CX3CR1', 'SPON2', 'GZMK','SPTSSB', 'IGFBP4','TNFRSF4', 'KIT', 'IL1R1','S100A8', 'S100A9', 'CD14','FCGR3A','FLT3', 'CD1C', 'CLEC10A',
# 'CLEC4C','C1QA', 'TREM2', 'APOE', 'CD9','MME', 'WASF1','CD79A', 'MS4A1', 'CD19','ZBTB32', 'ITGAX', 'TBX21','COCH', 'IGHA1', 'TNFRSF13B','TNFRSF17', 'MZB1', 'DERL3'
DefaultAssay(expData_s) = "RNA"
# datamat = expData_s$RNA@scale.data[which(rownames(expData_s) %in% markers_lst$gene), ]
datamat = expData_s$RNA@data[which(rownames(expData_s) %in% markers_lst$gene), ]
datamat = datamat[match(markers_lst$gene, rownames(datamat)), ]
datamat = datamat[, order(expData_s$cell_type)]
celltype = expData_s$cell_type
celltype = celltype[order(expData_s$cell_type)]

annotation_col = data.frame(cell_type = factor(celltype))      
rownames(annotation_col) = colnames(datamat)

# 1) rename cell_type 'column' to 'Cell type'
annotation_col <- annotation_col %>% 
  mutate('Cell type' = cell_type) %>% 
  select('Cell type')

# 1) clip heatmap values to +1
datamat[datamat >= 1.5] = 1.5

# 1) specify colors
mycolors <- c('#d60000',  '#ff7266', '#edb8b8', '#e36414',
  '#ff9f1c', '#ffbf69', '#f9cb9cff', '#fdf490', 
  '#a3b18a', '#588157', '#3a5a40', '#bce784',
  '#55a630', '#007f5f', '#a5be00', '#74d3ae',
  '#0774d8', '#00acc6')
# '#d60000',  '#ff7266', '#edb8b8', '#e36414',
# '#ff9f1c', '#ffbf69', '#f9cb9cff', '#fdf490', 
# '#a3b18a', '#588157', '#3a5a40', '#bce784',
# '#55a630', '#007f5f', '#a5be00', '#74d3ae',
# '#0774d8', '#00acc6'

# '#9651c4', '#bf03b8', '#eb0077', '#ff7ed1',  '#cdb4db', '#bcb6ff', 
# '#595959', '#a5a5a5', '#d2d4c8',  
# '#684756','#b69121', '#96705b', '#b76935', '#e7bc91'
names(mycolors) <- unique(annotation_col$'Cell type')
ann_colors <- list('Cell type' = mycolors)   

# 1) order cell_type
annotation_col$'Cell type'  = factor(annotation_col$'Cell type',
                                     levels = cell_type_list)

# 1) plot
mapal <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256)
datamat <- as.matrix(datamat)
p <- ComplexHeatmap::pheatmap(datamat,           # p <- 
              color = rev(mapal), 
              border_color = "grey60",
              name = "Expression",         # scale title
              annotation_col = annotation_col,
              annotation_colors = ann_colors,
              show_colnames = F, cluster_rows = F, cluster_cols = F,
              annotation_legend = F, 
              column_split = annotation_col$'Cell type', #gaps_col = T,
              column_title_gp = gpar(fontsize = 10), column_title_rot = 45,
              fontsize = 10, fontsize_row = 7)  

# 1) save
path_2 <- '/storage1/fs1/martyomov/Active/collaborations/romans/Greg_wu/CLOCK/output/figures/scRNA/UMAP_markers/'
save_pdf(p, paste0(path_2, "markers_T_mtx_v3.pdf"),  width = 5.5, height = 7)   #4.75


