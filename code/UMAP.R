#-------------------------------BCR/TCR: UMAP & CLONE SIZE ------------------------------

# 1) load files
path_5 <- '/storage1/fs1/martyomov/Active/collaborations/romans/Greg_wu/CLOCK/'
setwd(path_5)
expData <- readRDS('output/output/scRNA/combined_annot.rds')
contig_bcr_f <- read.csv('output/output/BCR/immc_H_subj_clons_germ_f.csv')
contig_tcr_f <- read.csv(paste0(path, "output/output/TCR/all_filtered_TCR.csv"))

## BCR mapping

# 1) count # unique clones 
uniq_clones <- contig_bcr_f$clone_id %>% 
  table() %>% as.data.frame() %>% 
  rename(clonotype_id = '.', clone_size = Freq)

# 1) remap clone_id col to expData
expData@meta.data <- expData@meta.data %>% select(-clonotype_id)
contig_annot_map <- contig_bcr_f %>% select(barcode, clone_id)
expData@meta.data <- merge(expData@meta.data, contig_annot_map, 
                           by=c("barcode"), all.x=TRUE)

# 1) map clone_size to expData
expData@meta.data <- expData@meta.data %>% rename(clonotype_id = clone_id)
expData@meta.data <- merge(expData@meta.data, uniq_clones, 
                           by=c("clonotype_id"), all.x=TRUE)
rownames(expData@meta.data) <- expData@meta.data$barcode
expData@meta.data <- expData@meta.data[order(expData@meta.data$ordered_id), ]
saveRDS(expData, 'output/output/scRNA/combined_annot.rds')

# 1) rename BCR-related colms
expData@meta.data <- expData@meta.data %>% rename(clone_size_bcr = clone_size,
                                                  clonotype_id_bcr = clonotype_id,
                                                  clone_group_bcr = clone_group)

## TCR mapping

# 1) count # unique clones 
uniq_clones <- contig_tcr_f$clonotype_id %>% 
  table() %>% as.data.frame() %>% 
  rename(clonotype_id_tcr = '.', clone_size_tcr = Freq)

# 1) remap clone_id col to expData
contig_annot_map <- contig_tcr_f %>% select(barcode, clonotype_id) %>% rename(clonotype_id_tcr = clonotype_id)
expData@meta.data <- merge(expData@meta.data, contig_annot_map, 
                           by=c("barcode"), all.x=TRUE)

# 1) map clone_size to expData
expData@meta.data <- merge(expData@meta.data, uniq_clones, 
                           by=c("clonotype_id_tcr"), all.x=TRUE)
rownames(expData@meta.data) <- expData@meta.data$barcode
expData@meta.data <- expData@meta.data[order(expData@meta.data$ordered_id), ]
saveRDS(expData, 'output/output/scRNA/combined_annot.rds')

# 1) plot "clone_size" using continius scale
expData@meta.data$clone_size_tcr[is.na(expData@meta.data$clone_size_tcr)] <- 0
FeaturePlot(expData, features = "clone_size_tcr")

# 1) add binned groups as I did with donuts plots
expData@meta.data <- expData@meta.data %>% 
  mutate(clone_group_tcr = case_when(clone_size_tcr == 0 ~ '0',
                                 clone_size_tcr == 1 ~ '1',
                                 clone_size_tcr >=2 & clone_size_tcr <= 4 ~ '2-4',
                                 clone_size_tcr >=5 & clone_size_tcr <= 9 ~ '5-9',
                                 clone_size_tcr >= 10 & clone_size_tcr <= 99 ~ '10-99',
                                 clone_size_tcr >= 100 ~ '>=100'))

# 1) fix the order & set the palette
expData@meta.data$clone_group_tcr = factor(expData@meta.data$clone_group_tcr,
                                   levels=c("0", "1", "2-4", "5-9", "10-99", ">=100"))
colors <- c("lightgray", "#ffccd5", "#ff8fa3", "#ff4d6d", "#d00000", "#6a040f")   #a4133c - 10

# 1) plot & save
DimPlot(expData, reduction = 'umap', pt.size = 0.15,
        group.by = 'clone_group_tcr', cols = colors) +
  labs(color = "Clone size") +
  ggtitle('T cell clonal expansion')

ggsave('output/figures/TCR/UMAP/UMAP_clone_size.png', width = 10, height = 7)
ggsave('output/figures/TCR/UMAP/UMAP_clone_size.pdf', width = 7, height = 5)



