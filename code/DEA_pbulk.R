library(reticulate)
use_python('/storage1/fs1/martyomov/Active/collaborations/romans/Greg_wu/CSF_atlas/code/scanpy/figs_de/r-reticulate/bin/python3.8')
library(Seurat)
library(anndata)
library(dplyr)
library(tidyverse)
library(viridis)
library(data.table)
library(ggrepel)

path <- "/storage1/fs1/martyomov/Active/collaborations/romans/Greg_wu/CLOCK/output/output/"
setwd(path)

# 1. upload RDS obj
seurat_data <- readRDS('scRNA/combined_annot.rds')
seurat_data_s <- seurat_data

# 2.1. filter out useless colms (+replace Proliferative T/NK)
org <- 'PBMC'
t <- 'YR2'
seurat_data@meta.data <- seurat_data_s@meta.data %>% 
  filter(tissue %in% org & time %in% c('BSE', t)) %>%   # 'W05', 'W10'
  mutate(time = recode(time, 'W05' = 'W5_10', 'W10' = 'W5_10')) %>%
  mutate(cell_type = recode(cell_type, 'Proliferative T/NK' = 'Proliferative T&NK'))

# 2.2. for CSF_W5_10 only:
seurat_data@meta.data <- seurat_data@meta.data %>% 
  filter(!cell_type %in% 'B Transitional')
# 2.3. for PBMC only:
seurat_data@meta.data <- seurat_data@meta.data %>% 
  filter(!cell_type %in% 'CD8 Trm')

# 3. create pseudobulk
pbulk <- AggregateExpression(
  object = seurat_data,assays = 'RNA',
  features = NULL,  return.seurat = TRUE,
  group.by = c('cell_type', 'sample'),
  add.ident = NULL, slot = "counts", verbose = TRUE)
pbulk@assays$RNA@counts %>% max()

# 4. normilize (for DE)
pbulk <- NormalizeData(pbulk, normalization.method = "LogNormalize", scale.factor = 10000)

# 5. scale (for heatmap)
all.genes <- rownames(pbulk)
pbulk <- ScaleData(pbulk, features = all.genes)

# 6. modify data-2: 2nd line is optional
pbulk@meta.data$cell_type.subject.time.organ <- rownames(pbulk@meta.data)
pbulk@meta.data <- pbulk@meta.data %>%
  mutate(time = sapply(strsplit(as.character(cell_type.subject.time.organ), split = "_"), `[`, 3)) %>%
  # mutate(time = recode(time, 'W05' = 'W5_10', 'W10' = 'W5_10')) %>%
  rename(cell_type = orig.ident) %>% 
  mutate(cell_type.time = paste0(cell_type, '_', time))

pbulk@meta.data$cell_type.time %>% unique()


# 7. reassign Idents to "cell_type.organ" (for DE: comparation)
Idents(pbulk) <- pbulk[['cell_type.time']]   
levels(pbulk) 

# 8.1. DE
for (i in 1:length(unique(pbulk$cell_type))){
  ident_1 <- paste0(unique(pbulk$cell_type)[i], "_BSE")
  ident_2 <- paste0(unique(pbulk$cell_type)[i], "_", t)
  pbulk.markers <- FindMarkers(pbulk, 
                               ident.1 = ident_1, ident.2 = ident_2,
                               test.use = "DESeq2", slot = "counts",
                               min.cells.group = 1, min.pct = 0.1,     # only.pos = T,
                               logfc.threshold = 0.1)
  pbulk.markers <- pbulk.markers %>% 
    mutate(cell_type = unique(pbulk$cell_type)[i])
  write.csv(pbulk.markers, file = paste0(path, '/scRNA/DEA/', org, '_', t, '/', unique(pbulk$cell_type)[i], '.csv'))}

# 8.2. merge DEGs from all cell types into 1 file
setwd(paste0(path, '/scRNA/DEA/', org, '_', t))
files <- list.files(pattern = ".csv")
combined_files <- bind_rows(lapply(files, fread))
write.csv(combined_files, file = paste0('combined.csv'))

# 9. load & modify data a bit
CSF_W5_10 <- read.csv('combined.csv')
CSF_W5_10 <- CSF_W5_10 %>% 
  rename(gene = V1) %>%
  subset(select = -X) %>% 
  mutate(log10_pval_log2FC = -log10(p_val)*avg_log2FC)

# CSF_W5_10 <- CSF_W5_10 %>% filter(cell_type %in% 'CD4 Naive')

# 10. add annotation to Color and Fill
CSF_W5_10 <- CSF_W5_10 %>% 
  mutate(gene_type = ifelse(p_val_adj < 0.05, "p.adj < 0.05", "p.adj > 0.05"))

Color <- c("p.adj < 0.05" = "#c1121f", "p.adj > 0.05" = "#495057")
Fill <- c("p.adj < 0.05" = "#c1121f", "p.adj > 0.05" = "white")

# 11. fix the order
CSF_W5_10$gene_type <- factor(CSF_W5_10$gene_type,
                              levels = c("p.adj < 0.05", "p.adj > 0.05"))

# 12. fix the order of labels
CSF_W5_10$cell_type = factor(CSF_W5_10$cell_type,
                         levels=c('CD4 Naive', 'CD4 Treg naive', 'CD4 Treg memory', 'CD4 Th1',
                                  'CD4 Th17', 'CD4 Th22', 'CD4 Temra', 'CD4 exhausted',
                                  'CD8 Naive', 'CD8 Tcm CCR4-', 'CD8 Tcm CCR4+', 
                                  #'CD8 Trm',
                                  'CD8 Tem GZMK+', 'CD8 Tem GZMB+', 'CD8 Temra', 'CD8 NKT-like',
                                  'gdT2', 'MAIT', 'Proliferative T&NK', 'NK CD56dim',
                                  'NK CD56bright', 'ILC', 'CD14 Mono', 'CD16 Mono',
                                  'cDC', 'pDC', 
                                  #'Microglia-like', 
                                  'B Transitional',
                                  'B Naive', 'B Atypical', 'B Switched', 'Plasmablast'))   # 'Microglia-like',

# 8. create a col "labels" 
CSF_W5_10$labels <- NA
CSF_W5_10$labels <- ifelse(CSF_W5_10$gene_type == 'p.adj > 0.05', NA, CSF_W5_10$gene)

# 13. plot
CSF_W5_10 %>% 
  ggplot(aes(x = cell_type, y = log10_pval_log2FC, 
             colour = gene_type, fill = gene_type)) +  #label = labels
  geom_jitter(shape=21, alpha = 0.5, size=1, width = 0.1) +
  # geom_text_repel(size = 2,  box.padding = 0.25,
  #                  max.overlaps = 50, aes(segment.size=0.2)) +
  geom_hline(yintercept = -0, color = "495057") +
  theme_classic() +
  scale_color_manual(values = Color) + 
  scale_fill_manual(values = Fill) +  
  xlab("Cell type") +
  ylab("-log10(p-value) * sign log2FC") +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  theme(legend.title=element_blank(),
        legend.position='top') + 
  guides(fill = guide_legend(nrow = 1)) 

# 14. save
path_2 <- "/storage1/fs1/martyomov/Active/collaborations/romans/Greg_wu/CLOCK/output/"
ggsave(paste0(path_2, "figures/scRNA/DEGs_manhattan_plot/", org, "_", t, ".png"),  #_lab
       width = 9, height = 5, dpi = 300)
# ggsave(paste0(path_2, "figures/scRNA/DEGs_manhattan_plot/", org, "_", t, "_lab.pdf"),
#        width = 9, height = 5, dpi = 300)
