library(reticulate)
use_python('/storage1/fs1/martyomov/Active/collaborations/romans/Greg_wu/CSF_atlas/code/scanpy/figs_de/r-reticulate/bin/python3.8')
library(Seurat)
library(anndata)
library(dplyr)
library(tidyverse)
library(viridis)

path <- "/storage1/fs1/martyomov/Active/collaborations/romans/Greg_wu/CLOCK/"
setwd(paste0(path, 'output/output/'))

# 1. upload raw h5ad obj
adata <- read_h5ad(paste0("combined_annot.h5ad"))
seurat_data <- CreateSeuratObject(counts = t(as.matrix(adata$layers['counts'])),      
                                  meta.data = adata$obs)

# 2. filter out useless colms
seurat_data@meta.data <- seurat_data@meta.data %>% 
  filter(tissue %in% 'PBMC' & time %in% c('BSE', 'W05', 'W10')) %>% 
  mutate(time = recode(time, 'W05' = 'W5_10', 'W10' = 'W5_10'))

# 3. create pseudobulk
pbulk <- AggregateExpression(
  object = seurat_data,assays = 'RNA',
  features = NULL,  return.seurat = TRUE,
  group.by = c('cell_type', 'sample'),
  add.ident = NULL, slot = "counts", verbose = TRUE)

pbulk@assays$RNA@counts %>% max()

# 5. normilize (for DE)
pbulk <- NormalizeData(pbulk, normalization.method = "LogNormalize", scale.factor = 10000)

# 6. scale (for heatmap)
all.genes <- rownames(pbulk)
pbulk <- ScaleData(pbulk, features = all.genes)

# 7. modify data-2: 2nd line is optional
pbulk@meta.data$cell_type.subject.time.organ <- rownames(pbulk@meta.data)
pbulk@meta.data <- pbulk@meta.data %>%
  mutate(time = sapply(strsplit(as.character(cell_type.subject.time.organ), split = "_"), `[`, 3)) %>%
  mutate(time = recode(time, 'W05' = 'W5_10', 'W10' = 'W5_10')) %>% 
  rename(cell_type = orig.ident) %>% 
  mutate(cell_type.time = paste0(cell_type, '_', time))

pbulk@meta.data$cell_type.time %>% unique()

# 8. reassign Idents to "cell_type.organ" (for DE: comparation)
Idents(pbulk) <- pbulk[['cell_type.time']]   
levels(pbulk) <- c("B Atypical_BSE","B Atypical_W5_10","B Naive_BSE","B Naive_W5_10"           
                   ,"B Switched_BSE","B Switched_W5_10","B Transitional_BSE","B Transitional_W5_10"    
                   ,"CD4 Naive_BSE","CD4 Naive_W5_10","CD4 Temra_BSE","CD4 Temra_W5_10"         
                   ,"CD4 Th1_BSE","CD4 Th1_W5_10","CD4 Th17_BSE","CD4 Th17_W5_10"          
                   ,"CD4 Th22_BSE","CD4 Th22_W5_10","CD4 Treg memory_BSE","CD4 Treg memory_W5_10"   
                   ,"CD4 Treg naive_BSE","CD4 Treg naive_W5_10","CD4 exhausted_BSE","CD4 exhausted_W5_10"     
                   ,"CD8 Naive_BSE","CD8 Naive_W5_10","CD8 Tcm CCR4-_BSE","CD8 Tcm CCR4-_W5_10"     
                   ,"CD8 Tem GZMB+_BSE","CD8 Tem GZMB+_W5_10","CD8 Tem GZMK+_BSE","CD8 Tem GZMK+_W5_10"     
                   ,"CD14 Mono_BSE","CD14 Mono_W5_10","CD16 Mono_BSE","CD16 Mono_W5_10"         
                   ,"ILC_BSE","ILC_W5_10","MAIT_BSE","MAIT_W5_10","NK CD56bright_BSE"       
                   ,"NK CD56bright_W5_10","NK CD56dim_BSE","NK CD56dim_W5_10","Plasmablast_BSE"         
                   ,"Plasmablast_W5_10","Proliferative T/NK_BSE","Proliferative T/NK_W5_10",
                   "cDC_BSE","cDC_W5_10","gdT1_BSE","gdT1_W5_10"              
                   ,"gdT2_BSE","gdT2_W5_10","pDC_BSE","pDC_W5_10")
levels(pbulk) 

# 9. DE
# for MAST only (Bug): rename cell names
cnames <- paste0('test_', rownames(pbulk@meta.data))
pbulk <- RenameCells(pbulk, new.names = cnames)
pbulk$cell_type.time %>% unique()

# cannot run pbulk because W5_10 (PBMC) and YR1_2 (CSF) groups
# will have clusters with less than 3 donors ("cells") per cluster
for (i in 1:length(unique(pbulk$cell_type))){
  ident_1 <- paste0(unique(pbulk$cell_type)[i], "_BSE")
  ident_2 <- paste0(unique(pbulk$cell_type)[i], "_W5_10")
  pbulk.markers <- FindMarkers(pbulk, 
                               ident.1 = ident_1,
                               ident.2 = ident_2,
                               test.use = "MAST", # MAST or wilcox?
                               slot = "data",
                               min.pct = 0.1,     # only.pos = T,
                               logfc.threshold = 0.1)
  write.csv(pbulk.markers, file = paste0('DEA/', unique(pbulk$cell_type)[i], '.csv'))
}



# 10. filter genes based on p_val_adj
pbulk.markers <- pbulk.markers %>% 
  rename(cell_type = cluster)   
# filter(p_val_adj < 0.05) 

# 11. visualize top_n markers per group/cell type
top_n <- pbulk.markers %>%
  group_by(cell_type) %>%
  top_n(n = 10, wt = avg_log2FC)

# top_n_2 <- pbulk.markers %>%
#   filter(gene %in% 'AXL')
# 
# top_n <- rbind(top_n_2, top_n)

levels(pbulk) <- c("CD32B+ cDC2_CSF", "CD32B+ cDC2_PBMC", 
                   "AREG+ cDC2_CSF",    "CD36+ cDC2_CSF",      
                   "CD36+ cDC2_PBMC")
pbulk@meta.data$cell_type.organ <- factor(pbulk@meta.data$cell_type.organ,
                                          levels = c("CD32B+ cDC2_CSF", "CD32B+ cDC2_PBMC", 
                                                     "AREG+ cDC2_CSF",    "CD36+ cDC2_CSF",      
                                                     "CD36+ cDC2_PBMC"))
DoHeatmap(pbulk, 
          features = top_n$gene,
          group.by = "cell_type.organ",
          size = 8,    
          angle = 45,
          draw.lines = T) +
  scale_fill_viridis(option = 'inferno') +                  
  guides(color = FALSE) +                   # hide cell_type legend
  theme(text = element_text(size = 16))

ggsave(paste0(path, "figures/myeloid/whole/new_DC2/CSF_PBMC_pbulk_DE_10_MAST.pdf"),
       width = 20,    #12
       height = 12,    #6
       dpi = 300)

