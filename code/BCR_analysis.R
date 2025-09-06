library(ggplot2)  # restart if R doesn't see ‘3.4.0’
library(dplyr)
library(stringr)
library(Seurat)
library(forcats)
library(mgsub)
library(DescTools)
library(ggpubr)
library(rstatix)
library(tidyverse)
library(ggbreak) 
# 1) load all immcantation specific packages
library(alakazam)
library(shazam)
library(tigger)
library(scoper)
library(dowser)
library(ggtree)
# library(plyr)

# -----------------------------COMBINE ALL BCRs INTO ONE DATAFRAME----------------

path_2 <- "/storage1/fs1/martyomov/Active/collaborations/romans/Greg_wu/CLOCK/data/BCR"
setwd(path_2)

merged_df <- list.files(pattern = 'filtered_contig_annotations\\.csv', recursive = TRUE) %>% 
  lapply(function(lst_file) {
    contig_df <- read.csv(lst_file)
    contig_df['patient_id'] <- strsplit(lst_file, split = '_', fixed=T)[[1]][2]
    contig_df['time'] <- strsplit(lst_file, split = '_', fixed=T)[[1]][3]
    contig_df['tissue'] <- strsplit(lst_file, split = '_', fixed=T)[[1]][4]
    contig_df['tissue'] <- ifelse(grepl('P', contig_df['tissue']), 'PBMC', 'CSF')
    contig_df['sample_id'] <- paste0(contig_df$patient_id, '_', 
                                     contig_df$time, '_', contig_df$tissue)
    contig_df['unique_barcode'] <- paste0(contig_df$sample_id, '_', contig_df$barcode)
    return(contig_df)
  }) %>%
  bind_rows

path_3 <- '/storage1/fs1/martyomov/Active/collaborations/romans/Greg_wu/CLOCK/output/output/BCR/'
write.csv(merged_df, paste0(path_3, 'all_unfiltered_BCR.csv'))

#----------------------------FILTERING OF BCRs & UMAP-------------------------

# 1) set env & load the data
path <- "/storage1/fs1/martyomov/Active/collaborations/romans/Greg_wu/CLOCK/"
setwd(path)
all_contigs <- read.csv(paste0(path, "output/output/BCR/all_unfiltered_BCR.csv"))
airr_rear <- read_tsv(paste0(path, "data/BCR/1307_0015_BSE_C/airr_rearrangement.tsv"))

# 1) filter out non-productive chains 
# 2) trim "-1" at the end of a barcode
# 3) remove columns that are not useful
all_contigs_filt <- all_contigs %>% 
  mutate(productive = recode(productive,  'true' = 'True')) %>% 
  filter(productive %in% 'True') %>% 
  mutate(barcode = gsub("-1", "", unique_barcode)) %>% 
  select(barcode, chain, cdr3)

# 1) create a column with aggregated alpha & beta chain by barcode
all_contigs_filt <- all_contigs_filt %>% 
  mutate(chain.cdr3_aa = paste0(chain, ':', cdr3)) %>% 
  select(barcode, chain.cdr3_aa) %>%  
  group_by(barcode) %>% 
  mutate(cdr3s_aa = paste0(chain.cdr3_aa, collapse = ';')) %>% 
  select(-chain.cdr3_aa) %>% 
  unique()                   # option 2: distinct(barcode, .keep_all = TRUE)

# 1) count #chains & leave only barcodes with 1 heavy: 1 light
all_contigs_filt$n_heavy <- sapply(1:nrow(all_contigs_filt), function(i) str_count(all_contigs_filt$cdr3s_aa[i], "IGH:"))
all_contigs_filt$n_light_k <- sapply(1:nrow(all_contigs_filt), function(i) str_count(all_contigs_filt$cdr3s_aa[i], "IGK:"))
all_contigs_filt$n_light_l <- sapply(1:nrow(all_contigs_filt), function(i) str_count(all_contigs_filt$cdr3s_aa[i], "IGL:"))
all_contigs_filt <- all_contigs_filt %>% filter(n_heavy == 1 & n_light_k == 1 | n_heavy == 1 & n_light_l == 1)

# 1) extract cell IDs from GEX
# 2) remove BCRs that are not in GEX
expData <- readRDS('output/output/scRNA/combined_annot.rds')
bcr_present <- rownames(expData@meta.data)
all_contigs_filt <- all_contigs_filt %>% 
  filter(barcode %in% bcr_present)

# 1) save the checkpoint & load it again
write.csv(all_contigs_filt, 'output/output/BCR/all_filtered_BCR.csv')
all_contigs_filt <- read.csv('output/output/BCR/all_filtered_BCR.csv')
expData <- readRDS('output/output/scRNA/combined_annot.rds')

# 1) add all columns from BCR data to Seurat obj
# 2) mark the row order to reverse it after merge function
expData@meta.data <- expData@meta.data %>%
  select(-barcode) %>% rename(barcode = unique_barcode)
all_contigs_filt <- all_contigs_filt %>% rename(cdr3s_aa_bcr = cdr3s_aa)
expData@meta.data$ordered_id  <- 1:nrow(expData@meta.data)
expData@meta.data <- merge(expData@meta.data, all_contigs_filt, 
                           by=c("barcode"), all.x=TRUE)
rownames(expData@meta.data) <- expData@meta.data$barcode
expData@meta.data <- expData@meta.data[order(expData@meta.data$ordered_id), ]

# 1) add BCR data to Seurat obj
# 2) filter out non-B cells which have BCRs & save GEX
# 3) do the same for BCR table
match.index = match(Cells(expData), all_contigs_filt$barcode)
expData@meta.data$contains_bcr <- !is.na(match.index)
expData@meta.data <- expData@meta.data %>%
  filter(!(cell_type %in% c('NK CD56dim','NK CD56bright', 'ILC', 'CD14 Mono', 'CD16 Mono', 'cDC', 'pDC', 'Microglia-like', 'gdT2', 'MAIT', 'CD4 Naive', 'CD8 Tcm CCR4-', 'CD4 Th1', 'CD4 Th17', 'CD4 Th22', 'CD4 Treg naive', 'CD8 Naive', 'CD8 Tcm CCR4+', 'CD8 NKT-like', 'CD8 Tem GZMK+', 'CD4 Treg memory', 'CD4 Temra', 'CD4 exhausted', 'CD8 Temra', 'CD8 Tem GZMB+', 'CD8 Trm', 'Proliferative T/NK')
           & contains_bcr == TRUE))
cellsToKeep <- rownames(expData@meta.data)
write.csv(cellsToKeep, 'output/output/BCR/cellsToKeep.csv')
expData <- expData[, colnames(expData) %in% cellsToKeep]
saveRDS(expData, 'output/output/scRNA/combined_annot.rds')
all_contigs_filt <- all_contigs_filt %>% filter(barcode %in% cellsToKeep)
write.csv(all_contigs_filt, 'output/output/BCR/all_filtered_BCR.csv')

# 1) visualize cells with TCRs on the UMAP
expData_semi <- readRDS('output/output/scRNA/semi_filtered/combined_annot.rds')
highlighted_cells <- Cells(expData)[which(expData$contains_bcr)]
Idents(expData) <- expData$cell_type
DimPlot(expData, reduction = 'umap', cells.highlight = highlighted_cells, label.size = 2.5,
        repel=TRUE, label = TRUE, cols = 'grey', pt.size = 0.1, sizes.highlight = 0.3) +
  scale_color_manual(labels = c("Non-clonal", "Clonal"),
                     values = c("grey", "#d90429")) 
ggsave('output/figures/BCR/UMAP/bcr_umap.png', width = 10, height = 7)
ggsave('output/figures/BCR/UMAP/bcr_umap.pdf', width = 7, height = 5)

#---------------------ISOTYPE FREQUENCY ANALYSIS--------------------------------

##--------------- BASIC: TOTAL BY CELL TYPE
path <- "/storage1/fs1/martyomov/Active/collaborations/romans/Greg_wu/CLOCK/"
setwd(path)
# all_contigs <- read.csv(paste0(path, "output/output/BCR/all_unfiltered_BCR.csv"))
all_contigs_filt <- read.csv('output/output/BCR/all_filtered_BCR.csv')
expData <- readRDS('output/output/scRNA/combined_annot.rds')

# 1) get "c_gene" column together with barcodes
# 2) exclude cells without annotated isotype
all_filt_heavy <- all_contigs %>% 
  mutate(productive = recode(productive,  'true' = 'True')) %>% 
  filter(productive %in% 'True') %>% 
  mutate(barcode = gsub("-1", "", unique_barcode)) %>% 
  filter(barcode %in% all_contigs_filt$barcode & chain %in% "IGH") %>% 
  select(barcode, c_gene) %>% 
  rename(isotype = c_gene) %>% 
  filter(!isotype %in% "")

# 1) exclude the cell without annotated isotype in GEX
# 2) and then in BCR table
expData@meta.data <- expData@meta.data %>%
  filter(!barcode %in% "0015_BSE_PBMC_TCGCGTTGTGAAAGAG")
cellsToKeep <- rownames(expData@meta.data)
write.csv(cellsToKeep, 'output/output/BCR/cellsToKeep.csv')
expData <- expData[, colnames(expData) %in% cellsToKeep]
saveRDS(expData, 'output/output/scRNA/combined_annot.rds')
all_contigs_filt <- all_contigs_filt %>% filter(barcode %in% cellsToKeep)
write.csv(all_contigs_filt, 'output/output/BCR/all_filtered_BCR.csv')

# 1) add "isotype" column from BCR data to Seurat obj
# 2) mark the row order to reverse it after merge function
expData@meta.data <- merge(expData@meta.data, all_filt_heavy, 
                           by=c("barcode"), all.x=TRUE)
rownames(expData@meta.data) <- expData@meta.data$barcode
expData@meta.data <- expData@meta.data[order(expData@meta.data$ordered_id), ]

# 1) select cells with isotype in GEX
isotype_df <- expData@meta.data %>% 
  filter(contains_bcr == TRUE) %>% 
  select(cell_type, isotype)            # subject, time, tissue, sample,

# 1) calculate % of isotype usage for each cell_type 
isotype_df <- isotype_df %>%
  group_by(cell_type, isotype) %>% 
  summarise(per_cluster = n()) %>% 
  mutate(sum_cells = sum(per_cluster)) %>% 
  mutate(percent = per_cluster/sum_cells * 100)

# 1) CHECK
isotype_df %>% group_by(cell_type) %>% 
  summarise(sum_perc = sum(percent))

# 1) fix the order of labels 
isotype_df$cell_type = factor(isotype_df$cell_type,
                              levels=c('B Transitional', 'B Naive', 'B Atypical', 'B Switched', 'Plasmablast'))
isotype_df$isotype = factor(isotype_df$isotype,
                            levels=c("IGHM", "IGHD", "IGHA1", "IGHA2",
                                     "IGHG1", "IGHG2", "IGHG3", "IGHG4"))

# 1) add custom palette
my_palette <- c('#8ecae6', '#118ab2', '#ffc971', '#ff9505',
                '#ffa69e', '#ff4d6d', '#c9184a', '#590d22')  #023e8a

# 1) plot
ggplot(isotype_df, aes(x = cell_type, y = percent, 
                       fill = as.factor(isotype), label = isotype))+    
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = my_palette, name = 'Isotype') +   # uncomment & add name = 'Cell type' if u want to use custom palette
  ylab("% of cell barcodes") +
  theme_light() +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        legend.text=element_text(size=12),
        legend.title=element_text(size=13)) +
  scale_x_discrete(guide = guide_axis(angle = 45))

ggsave('output/figures/BCR/Isotype_freq/all_cell_types.png', width = 7, height = 5)
ggsave('output/figures/BCR/Isotype_freq/all_cell_types.pdf', width = 7, height = 5)

#--------------- ADVANCED: BY CELL TYPE & ORGAN/TIME----------------------------

## ORGAN
# 1) select cells with isotype in GEX
isotype_df <- expData@meta.data %>% 
  filter(contains_bcr == TRUE & time %in% "BSE") %>% 
  select(cell_type, tissue, isotype)            # subject, time, tissue, sample,

# 1) calculate % of isotype usage for each cell_type 
isotype_df <- isotype_df %>%
  group_by(cell_type, tissue, isotype) %>% 
  summarise(per_cluster = n()) %>% 
  mutate(sum_cells = sum(per_cluster)) %>% 
  mutate(percent = per_cluster/sum_cells * 100)

# 1) CHECK
isotype_df %>% group_by(cell_type, tissue) %>% 
  summarise(sum_perc = sum(percent))

# 1) fix the order of labels 
isotype_df$cell_type = factor(isotype_df$cell_type,
                              levels=c('B Transitional', 'B Naive', 'B Atypical', 'B Switched', 'Plasmablast'))
isotype_df$tissue = factor(isotype_df$tissue, levels=c('PBMC', 'CSF'))
isotype_df$isotype = factor(isotype_df$isotype,
                            levels=c("IGHM", "IGHD", "IGHA1", "IGHA2",
                                     "IGHG1", "IGHG2", "IGHG3", "IGHG4"))

# 1) add custom palette
my_palette <- c('#8ecae6', '#118ab2', '#ffc971', '#ff9505',
                '#ffa69e', '#ff4d6d', '#c9184a', '#590d22')  #023e8a

# 1) plot using stacked barplot
ggplot(isotype_df) +
  geom_bar(aes(x = tissue, y = percent, fill = isotype),
           position = "stack", stat = "identity", color = "black") +
  facet_grid(~cell_type) +          
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        axis.title.x = element_blank()) +
  ggtitle("Isotype frequencies in BSE: PBMC vs CSF") +
  ylab("% of cell barcodes") +
  scale_fill_discrete(name = "Isotype") +
  scale_fill_manual(values = my_palette, name = 'Isotype')  

ggsave('output/figures/BCR/Isotype_freq/BSE_PBMC_CSF_v2.png', width = 14, height = 5)
ggsave('output/figures/BCR/Isotype_freq/BSE_PBMC_CSF_v2.pdf', width = 14, height = 5)

## TIME
# 1) select cells with isotype in GEX
organ <- 'PBMC'
isotype_df <- expData@meta.data %>% 
  filter(contains_bcr == TRUE & tissue %in% organ) %>% 
  mutate(time = recode(time, 'W05' = 'W5_10', 'W10' = 'W5_10')) %>% 
  # filter(!time %in% "YR2") %>%             # for CSF only
  select(cell_type, time, isotype)            # subject, time, tissue, sample,

# 1) calculate % of isotype usage for each cell_type 
isotype_df <- isotype_df %>%
  group_by(cell_type, time, isotype) %>% 
  summarise(per_cluster = n()) %>% 
  mutate(sum_cells = sum(per_cluster)) %>% 
  mutate(percent = per_cluster/sum_cells * 100)

# 1) CHECK
isotype_df %>% group_by(cell_type, time) %>% 
  summarise(sum_perc = sum(percent))

# 1) fix the order of labels 
isotype_df$cell_type = factor(isotype_df$cell_type,
                              levels=c('B Transitional', 'B Naive', 'B Atypical', 'B Switched', 'Plasmablast'))
isotype_df$time = factor(isotype_df$time, levels=c('BSE', 'W5_10', 'YR1', 'YR2'))  #, 'YR2'
isotype_df$isotype = factor(isotype_df$isotype,
                            levels=c("IGHM", "IGHD", "IGHA1", "IGHA2",
                                     "IGHG1", "IGHG2", "IGHG3", "IGHG4"))

# 1) add custom palette
my_palette <- c('#8ecae6', '#118ab2', '#ffc971', '#ff9505',
                '#ffa69e', '#ff4d6d', '#c9184a', '#590d22')  
# for CSF only
my_palette <- c('#8ecae6', '#118ab2', '#ffc971',
                '#ffa69e', '#ff4d6d', '#c9184a', '#590d22')

# 1) plot using stacked barplot
ggplot(isotype_df) +
  geom_bar(aes(x = time, y = percent, fill = isotype),
           position = "stack", stat = "identity", color = "black") +
  facet_grid(~cell_type) +          
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        axis.title.x = element_blank()) +
  ggtitle(paste0("Isotype frequencies in ", organ)) +
  ylab("% of cell barcodes") +
  scale_fill_discrete(name = "Isotype") +
  scale_fill_manual(values = my_palette, name = 'Isotype')  

ggsave(paste0('output/figures/BCR/Isotype_freq/', organ, '_time.png'), width = 11, height = 4)
ggsave(paste0('output/figures/BCR/Isotype_freq/', organ, '_time.pdf'), width = 11, height = 4)

#-----------------------------IMMCANT: MODIFY IMMCANTATIONS' FILES----------------------

path_2 <- "/storage1/fs1/martyomov/Active/collaborations/romans/Greg_wu/CLOCK/data/BCR/immcantation/changeo_10x/"
setwd(path_2)

# 1) add metadata columns and rename barcodes
list.files(pattern = 'filtered_contig_heavy_germ-pass\\.tsv', recursive = TRUE) %>% 
  lapply(function(lst_file) {
    contig_df <- read.csv(lst_file, sep = "\t")
    contig_df['patient_id'] <- strsplit(lst_file, split = '_', fixed=T)[[1]][2]
    contig_df['time'] <- strsplit(lst_file, split = '_', fixed=T)[[1]][3]
    contig_df['tissue'] <- strsplit(lst_file, split = '_', fixed=T)[[1]][4]
    contig_df['tissue'] <- ifelse(grepl('P', contig_df['tissue']), 'PBMC', 'CSF')
    contig_df['sample_id'] <- paste0(contig_df$patient_id, '_', 
                                     contig_df$time, '_', contig_df$tissue)
    contig_df['unique_barcode'] <- paste0(contig_df$sample_id, '_', contig_df$cell_id) 
    contig_df <- contig_df %>% mutate(barcode = gsub("-1", "", unique_barcode))
    dir.create(paste0('Redef_Clones/subj_', unique(contig_df['patient_id'])), showWarnings = FALSE)
    write.table(contig_df, paste0('Redef_Clones/subj_', unique(contig_df['patient_id']), '/heavy_germ_',
                                  unique(contig_df['sample_id']) , '.tsv'),
                row.names =F, quote=FALSE, sep='\t')
    return(contig_df)
  }) 

# 1) join files of *heavy-germ-subj-X.tsv within each subject/folder
path_3 <- "/storage1/fs1/martyomov/Active/collaborations/romans/Greg_wu/CLOCK/data/BCR/immcantation/DefineClones"
setwd(path_3)

list.dirs('.', recursive=FALSE) %>% 
  lapply(function(lst_folder) {
    setwd(path_3)
    setwd(lst_folder)
    merged <- ldply(list.files(), read_tsv)
    write.table(merged, paste0("merged.tsv"), row.names =F,
                quote=FALSE, sep='\t')
  }) 

# 1) merge all *heavy-germ-subj-X.tsv into one
merged_df <- list.files(pattern = 'merged_p_clone-pass\\.tsv', recursive = TRUE) %>% 
  lapply(function(lst_file) {
    contig_df <- read.csv(lst_file, sep = "\t")
    contig_df['clone_id'] <- paste0('subj_', contig_df$patient_id, '_clone_', contig_df$clone_id)
    return(contig_df)
  }) %>%
  bind_rows

path_4 <- '/storage1/fs1/martyomov/Active/collaborations/romans/Greg_wu/CLOCK/output/output/BCR/'
write.csv(merged_df, paste0(path_4, 'immcant_unfilt_clons.csv'))

# -----------------------------REMAP IMMCANTATION BCR CLONES & ISOTYPES TO GEX----------------------

# 1) load files
path_5 <- '/storage1/fs1/martyomov/Active/collaborations/romans/Greg_wu/CLOCK/'
setwd(path_5)
all_contigs_filt <- read.csv('output/output/BCR/all_filtered_BCR.csv')
immc_unfilt_clons <- read.csv('output/output/BCR/immcant_unfilt_clons.csv')
expData <- readRDS('output/output/scRNA/combined_annot.rds')

# 1) filter out cells which are present in all_filtered_BCR df 
immc_filt_clons <- immc_unfilt_clons %>% 
  filter(barcode %in% all_contigs_filt$barcode) %>% 
  rename(clonotype_id = clone_id)
write.csv(immc_filt_clons, 'output/output/BCR/immcant_filt_clons.csv')

# 1) remove useless/outdated columns
expData@meta.data <- expData@meta.data %>% 
  select(-c(X, cdr3s_aa_bcr, n_heavy, n_light_k, n_light_l, contains_bcr, isotype, leiden_3.0,       
            leiden_3.5, leiden_4.0, leiden_4.5, leiden_5.0, leiden_5.5,))

# 1) leave only "barcode", "clone_id" & "c_call" column in immc_filt_clons
# 2) map "clone_id" and "c_call" column to GEX
# 3) mark the row order to reverse it after merge function
immc_clons_s <- immc_filt_clons %>% select(barcode, clone_id, c_call) %>% rename(isotype = c_call)
expData@meta.data <- merge(expData@meta.data, immc_clons_s, 
                           by=c("barcode"), all.x=TRUE)
rownames(expData@meta.data) <- expData@meta.data$barcode
expData@meta.data <- expData@meta.data[order(expData@meta.data$ordered_id), ]

# 1) add "contains_bcr" & save GEX
match.index = match(Cells(expData), immc_clons_s$barcode)
expData@meta.data$contains_bcr <- !is.na(match.index)
expData@meta.data <- expData@meta.data %>% rename(clonotype_id = clone_id)
saveRDS(expData, 'output/output/scRNA/combined_annot.rds')

#-------------------------------UMAP & CLONE SIZE ------------------------------

# 1) load files
path_5 <- '/storage1/fs1/martyomov/Active/collaborations/romans/Greg_wu/CLOCK/'
setwd(path_5)
expData <- readRDS('output/output/scRNA/combined_annot.rds')
expData@meta.data %>% colnames()

# 1) count # unique clones 
contig_annot_short <- read.csv('output/output/BCR/immc_H_subj_clons_germ_f.csv')
uniq_clones <- contig_annot_short$clone_id %>% 
  table() %>% as.data.frame() %>% 
  rename(clonotype_id = '.', clone_size_tcr = Freq)

# 1) remap clone_id col to expData
expData@meta.data <- expData@meta.data %>% select(-clonotype_id)
contig_annot_map <- contig_annot_short %>% select(barcode, clone_id)
expData@meta.data <- merge(expData@meta.data, contig_annot_map, 
                           by=c("barcode"), all.x=TRUE)

# 1) map clone_size to expData
expData@meta.data <- expData@meta.data %>% rename(clonotype_id = clone_id)
expData@meta.data <- merge(expData@meta.data, uniq_clones, 
                           by=c("clonotype_id"), all.x=TRUE)
rownames(expData@meta.data) <- expData@meta.data$barcode
expData@meta.data <- expData@meta.data[order(expData@meta.data$ordered_id), ]
saveRDS(expData, 'output/output/scRNA/combined_annot.rds')

# 1) plot "clone_size" using continius scale
expData@meta.data$clone_size[is.na(expData@meta.data$clone_size)] <- 0
FeaturePlot(expData, features = "clone_size")

# 1) add binned groups as I did with donuts plots
expData@meta.data <- expData@meta.data %>% 
  mutate(clone_group = case_when(clone_size == 0 ~ '0',
                                 clone_size == 1 ~ '1',
                                 clone_size >=2 & clone_size <= 4 ~ '2-4',
                                 clone_size >=5 & clone_size <= 9 ~ '5-9',
                                 clone_size >= 10 ~ '>=10'))

# 1) fix the order & set the palette
expData@meta.data$clone_group = factor(expData@meta.data$clone_group,
                                       levels=c("0", "1", "2-4", "5-9", ">=10"))
colors <- c("lightgray", "#ffccd5", "#ff8fa3", "#ff4d6d", "#a4133c") 

# 1) select "highlighted_cells"
highlighted_cells <- Cells(expData)[which(expData$contains_bcr)]
Idents(expData) <- expData$cell_type

# 1) plot & save
DimPlot(expData, reduction = 'umap', pt.size = 0.15,
        group.by = 'clone_group', cols = colors) +
  labs(color = "Clone size") +
  ggtitle('B cell and Plasmablast clonal expansion')

ggsave('output/figures/BCR/UMAP/UMAP_clone_size.png', width = 10, height = 7)
ggsave('output/figures/BCR/UMAP/UMAP_clone_size.pdf', width = 7, height = 5)

#-----------------------CLONAL DYNAMICS USING SCREP-----------------------------

# 1) load the data
path <- "/storage1/fs1/martyomov/Active/collaborations/romans/Greg_wu/CLOCK/"
setwd(path)
all_contigs_filt <- read.csv('output/output/BCR/immc_H_subj_clons_germ_f.csv')
all_contigs_unf <- read.csv('output/output/BCR/all_unfiltered_BCR.csv')

# 1) modify & filter cells by barcode
contigs_full_f <- all_contigs_unf %>% 
  mutate(barcode = gsub("-1", "", unique_barcode)) %>% 
  filter(barcode %in% all_contigs_filt$barcode)

# 1) reformat data into scRep format of lists
path_2 <- "/storage1/fs1/martyomov/Active/collaborations/romans/Greg_wu/CLOCK/data/BCR"
setwd(path_2)

merged_df <- list.files(pattern = 'filtered_contig_annotations\\.csv', recursive = TRUE) %>% 
  lapply(function(lst_file) {
    contig_df <- read.csv(lst_file)
    contig_df['patient_id'] <- strsplit(lst_file, split = '_', fixed=T)[[1]][2]
    contig_df['time'] <- strsplit(lst_file, split = '_', fixed=T)[[1]][3]
    contig_df['tissue'] <- strsplit(lst_file, split = '_', fixed=T)[[1]][4]
    contig_df['tissue'] <- ifelse(grepl('P', contig_df['tissue']), 'PBMC', 'CSF')
    contig_df['sample'] <- paste0(contig_df$patient_id, '_', 
                                  contig_df$time, '_', contig_df$tissue)
    contig_df['unique_barcode'] <- paste0(contig_df$sample, '_', contig_df$barcode)
    contig_df <- contig_df %>% 
      mutate(unique_barcode = gsub("-1", "", unique_barcode)) %>% 
      filter(unique_barcode %in% all_contigs_filt$barcode) %>% 
      select(-unique_barcode)
    return(contig_df)
  })

path_3 <- '/storage1/fs1/martyomov/Active/collaborations/romans/Greg_wu/CLOCK/output/output/BCR/'
saveRDS(merged_df, paste0(path_3, "full_filtered_BCR_list.RData"))
BCR.contigs <- readRDS(paste0(path_3, 'full_filtered_BCR_list.RData'))


for (i in 1:length(BCR.contigs)){
  print(unique(full_filtered_TCR[[i]]$sample)) 
}

# 1) combine contigs into clones
combined.BCR <- combineBCR(BCR.contigs, 
                           samples = c("001_YR2_CSF","001_YR2_PBMC" ,"0015_BSE_CSF","0015_BSE_PBMC"
                                       ,"003_YR2_PBMC","004_YR2_PBMC" ,"005_BSE_CSF","005_BSE_PBMC"
                                       ,"005_YR1_CSF","005_YR1_PBMC" ,"006_BSE_CSF","006_BSE_PBMC"
                                       ,"006_YR1_CSF","006_YR1_PBMC" ,"006_YR2_PBMC","007_W10_PBMC"
                                       ,"009_BSE_PBMC","009_YR2_PBMC" ,"010_BSE_PBMC","010_YR2_PBMC"
                                       ,"011_BSE_CSF","011_W10_CSF" ,"011_W10_PBMC","011_YR2_PBMC"
                                       ,"012_BSE_PBMC","012_YR2_PBMC" ,"013_BSE_CSF","013_BSE_PBMC"
                                       ,"013_W05_CSF","013_W05_PBMC"),
                           threshold = 0.85,
                           removeNA = TRUE, 
                           removeMulti = TRUE, 
                           filterMulti = FALSE)

head(combined.BCR[[1]])

combined.BCR.df <- all_contigs_filt %>% 
  filter(sample_id %in% c("013_BSE_PBMC","013_W05_PBMC", "013_BSE_CSF", "013_W05_CSF"))

combined.BCR.df$sample_id = factor(combined.BCR.df$sample_id,
                               levels=c("013_BSE_PBMC","013_W05_PBMC", "013_BSE_CSF", "013_W05_CSF"))

# 1) plot alluvial stacked barplot
# "005_BSE_PBMC",  "005_YR1_PBMC", "005_BSE_CSF", "005_YR1_CSF"
# "006_BSE_PBMC", "006_YR1_PBMC", "006_YR2_PBMC", "006_BSE_CSF", "006_YR1_CSF"
# "011_BSE_PBMC" ,"011_W10_PBMC" ,"011_YR2_PBMC", "011_BSE_CSF" ,"011_W10_CSF"
# "013_BSE_PBMC","013_W05_PBMC", "013_BSE_CSF", "013_W05_CSF" 
p <- clonalCompare(combined.BCR, top.clones = 5, 
                   samples = c("013_BSE_PBMC", "013_BSE_CSF", "013_W05_CSF"), 
                   cloneCall="CTstrict", graph = "alluvial") 

p + theme(axis.text.x = element_text(angle=45, hjust=1)) +
  ylab("Proportion of clonotypes") +
  ggtitle("Top clonotypes in subj_13") +
  scale_x_discrete(limits = c("013_BSE_PBMC", "013_BSE_CSF", "013_W05_CSF"))
  

ggsave(paste0(path, 'output/figures/BCR/Clone_barplot/subj_13_PBMC_CSF.png'), 
       width = 5, height = 4)
ggsave(paste0(path, 'output/figures/BCR/Clone_barplot/subj_13_PBMC_CSF.pdf'), 
       width = 5, height = 4)

#---------------------- DOTPLOT WITH CLONOTYPES ACROSS ALL B CELLS--------------

# 1) count # unique clones 
contig_annot_short <- read.csv('output/output/BCR/immcant_heavy_clons_f.csv')
uniq_clones <- contig_annot_short$clonotype_id %>% 
  table() %>% as.data.frame()

# 1) compress it, by calculating clone size & rename columns
uniq_clones <- uniq_clones$Freq %>%
  table() %>% as.data.frame() %>% 
  rename(clone_size = '.', n_clonotypes = Freq)

# 1) plot 
ggplot(uniq_clones, aes(x = reorder(clone_size, -n_clonotypes), y = n_clonotypes)) +
  geom_bar(stat="identity") +
  theme_light() +
  xlab("Clone size") +
  ylab("# Clonotypes") +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  scale_y_break(c(69, 4650))     
ggsave('output/figures/BCR/Freq_barplot/n_clonotypes.png', width = 10, height = 7)

#-----------------------------------CLONAL DIVERSITY: DONUT CHART----------------------------

# 1) DONUT CHART: calculate frequncies of clonotypes grouped by tissue & cell_type_2 & pick the biggest
freq_bcr <- expData@meta.data %>%
  # filter(subject %in% c("006", "009", "010", "011", "012")) %>%   #filter matched samples
  mutate(time_comb = time) %>% 
  mutate(time_comb = recode(time_comb, 'W05' = 'W5_10', 'W10' = 'W5_10')) %>%   #'YR1' = 'YR1_2', 'YR2' = 'YR1_2',
  filter(contains_bcr == TRUE) %>% 
  filter(cell_type %in% c('B Naive', 'B Transitional', 'B Atypical', 'B Switched', 'Plasmablast')) %>%
  mutate(cell_type_2 = 'B cell') %>% 
  group_by(cell_type_2, tissue, time_comb, clonotype_id) %>% 
  tally() 

# 1) keep only those groups in which you are interested in
freq_bcr_filt <- freq_bcr %>% 
  filter(tissue %in% c('CSF') & 
           time_comb %in% c('BSE', 'YR1')) %>%   #, 'YR1_2'
  rename(clon_size = "n")
freq_bcr_filt$clon_size %>% table()

# 1) count the total #clones per time group
# 2) calculate in % the contribution of clones in each group
num_tot_cells <- freq_bcr_filt %>% 
  group_by(time_comb) %>%   #tissue; time_comb
  summarise(total_clon_size = sum(clon_size)) 
freq_bcr_filt$total_clon_size <- with(freq_bcr_filt, 
                                      num_tot_cells$total_clon_size[match(time_comb, num_tot_cells$time_comb)])   #time_comb; tissue
freq_bcr_filt <- freq_bcr_filt %>% 
  mutate(percent = clon_size/total_clon_size * 100)

# 1) create the colomn "clone_group" using when_case & split into 4 groups (unique, 2-4, 5-9, >=10)
freq_bcr_filt <- freq_bcr_filt %>% 
  mutate(clone_group = case_when(clon_size == 1 ~ 'Unique',
                                 clon_size >=2 & clon_size <= 4 ~ '2-4',
                                 clon_size >=5 & clon_size <= 9 ~ '5-9',
                                 clon_size >= 10 ~ '>=10'))

# 1) collapse based on clone_group/Unique group
freq_bcr_filt <- freq_bcr_filt %>% 
  group_by(cell_type_2, tissue, time_comb, 
           grp = case_when(clone_group %in% "Unique" ~ 0, TRUE ~ row_number())) %>% 
  summarise(across(starts_with(c("clon_size", "percent")), sum),
            clone_group = str_c(clone_group, collapse=",")) %>% 
  select(-grp) %>% 
  mutate(clone_group = gsub("^Uni.*", "Unique", clone_group))

# 1) order the labels
freq_bcr_filt$clone_group %>% unique()
freq_bcr_filt$clone_group = factor(freq_bcr_filt$clone_group,
                                   levels=c("Unique", "2-4", "5-9", ">=10"))
# freq_bcr_filt$tissue = factor(freq_bcr_filt$tissue,
#                                    levels=c("PBMC", "CSF"))

# plot the Donut chart
colors <- c("#adb5bd", "#ff8fa3", "#ff4d6d", "#a4133c")
ggplot(data = freq_bcr_filt, aes(x=2, y = percent))+    #, fill = clone_group
  geom_bar(stat="identity", color='black')+
  geom_bar(stat="identity", aes(fill = clone_group))+
  # geom_col(color='black')+
  coord_polar("y", start = 0) + 
  facet_wrap(~ time_comb) +    #tissue; time_comb
  theme(panel.background = element_blank(),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(), 
        strip.text = element_text(size = 15)) +
  scale_fill_manual(values = colors, name = 'No. of clones') +
  xlim(0.5, 2.5)  

ggsave('output/figures/BCR/Donut_chart/across_all/CSF_BSE_YR1_sep.png', width = 10, height = 7)
ggsave('output/figures/BCR/Donut_chart/across_all/CSF_BSE_YR1_sep.pdf', width = 10, height = 7)

# 1) count # B cells with BCR per tissue per time for the presa
n_bcr_tissue_time <- expData@meta.data %>% 
  mutate(time_comb = time) %>% 
  mutate(time_comb = recode(time_comb, 'W05' = 'W5_10', 'W10' = 'W5_10')) %>% 
  filter(contains_bcr == TRUE & tissue %in% "CSF") %>% 
  group_by(tissue, time_comb) %>% 
  tally()

#------------------------CLONAL DIVERSITY: ADVANCED BOXPLOT: BY TIMEPOINTS--------------------------------

org <- 'PBMC'
time_filt <- 'W5_10'
clons_for_gini <- expData@meta.data %>% 
  mutate(time = recode(time, 'W05' = 'W5_10', 'W10' = 'W5_10')) %>%   #'YR1' = 'YR1_2', 'YR2' = 'YR1_2',
  filter(tissue %in% org &             ##
           cell_type %in% c('B Naive', 'B Transitional', 'B Atypical', 'B Switched', 'Plasmablast') &
           time %in% c("BSE", time_filt)                        ## , time_filt
  ) %>%   ## subject %in% c("005", "006", "009", "011", "013", "0015")
  filter(contains_bcr == TRUE) %>% 
  group_by(cell_type, time, sample, clonotype_id_bcr) %>%    # tissue, time
  tally() %>% 
  rename(freq = "n")

# 1) calculate the Gini coefficient per cell_type & per sample
gini_res <- clons_for_gini %>% 
  group_by(cell_type, time, sample) %>%        ## tissue, time,
  summarise(across(-clonotype_id_bcr, Gini, na.rm = TRUE)) %>% 
  filter(!freq %in% NaN)

# 1) order the labels
# 2) specify the consistent colors
gini_res$cell_type = factor(gini_res$cell_type, levels=c('B Naive', 'B Transitional', 'B Atypical', 'B Switched', 'Plasmablast'))
gini_res$time = factor(gini_res$time, levels=c('BSE', time_filt))    
# gini_res$tissue = factor(gini_res$tissue, levels=c('PBMC', 'CSF'))

# 1) create a column "subject"
gini_res$sample <- as.character(gini_res$sample)
gini_res['subject'] <- sapply(gini_res$sample, function(x){
  strsplit(x, split = "_", fixed = T)[[1]][1]})

gini_res <- gini_res %>% 
  filter(!cell_type %in% c('B Transitional'))
# 1) exclude CD8 Trm because they are only in CSF
# gini_res <- gini_res %>% filter(!cell_type %in% "CD8 Trm")
# gini_res <- gini_res %>% filter(!cell_type %in% "CD8 NKT-like")   # for CSF: W5_10

# 1) run tukey test (NORMAL)
tukey_stat <- gini_res %>%
  group_by(cell_type) %>%
  tukey_hsd(freq ~ time, p.adjust.method = 'fdr') %>%     # tissue; time
  filter(group1 == 'BSE') %>%          # PBMC; BSE
  add_xy_position(x = "cell_type") 

# 1) specify your pallete
# BSE: #ea698b, W5_10: #dec9e9, YR1: #a06cd5, YR2: #47126b
my_palette <- c('#ea698b', '#a06cd5')    

# 1) BOXPLOTS WITH STATS: plot the values of Gini coefficient per cell type
fig <- ggplot(gini_res, aes(x = cell_type, y = freq)) +
  geom_boxplot(outlier.shape = NA, aes(fill = time), alpha = 0.7) +    # tissue; alpha = 0.5
  geom_point(position=position_dodge(width=0.75), alpha = 0.7,        # alpha = 0.5
             aes(fill = time), colour="black", pch=21) +             # tissue
  geom_line(aes(x = as.numeric(factor(cell_type)) + .75 / 4 * ifelse(time == "BSE", -1, 1),   # tissue == "PBMC"; time == "BSE"
                group = interaction(subject, cell_type)), color = "black", linewidth = 0.5, alpha = 0.5) +
  ylab("Gini coefficient") +
  theme_light() +
  theme(axis.title.x = element_blank(),
        axis.title.y=element_text(size=10),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        legend.position="right",
        legend.title=element_text(size=9),
        legend.text=element_text(size=8),
        plot.title = element_text(size = 10, hjust = 0.5))   +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  stat_pvalue_manual(tukey_stat, hide.ns = T) +
  # scale_fill_discrete(name = "Tissue") +
  scale_fill_manual(name = "Time", values = my_palette) +
  ggtitle(paste0("BCR repertoire diversity in ", org, ": BSE vs ", time_filt))   #"TCR repertoire diversity in BSE: PBMC vs CSF"

# 1) save in this to keep the good asterisk shape
ggsave(filename=paste0('output/figures/BCR/Gini_coef/total/', org, '_BSE_', time_filt, '.pdf'), 
       width = 5, height = 3, dpi = 300)
png(filename=paste0('output/figures/BCR/Gini_coef/total/', org, '_BSE_', time_filt, '.png'), 
    width = 5, height = 3, units = "in", res = 300)   # , org, '_BSE_', time_filt, '.png'
plot(fig)
dev.off()

#------------------------CLONAL DIVERSITY: ADVANCED BOXPLOT: BY TISSUE--------------------------------

org <- 'CSF'
time_filt <- 'YR1'
clons_for_gini <- expData@meta.data %>% 
  mutate(time = recode(time, 'W05' = 'W5_10', 'W10' = 'W5_10')) %>%   #'YR1' = 'YR1_2', 'YR2' = 'YR1_2',
  filter(#tissue %in% org &             ##
    cell_type %in% c('B Naive', 'B Transitional', 'B Atypical', 'B Switched', 'Plasmablast') &
      time %in% c("BSE")                        ## , time_filt
  ) %>%   ## subject %in% c("005", "006", "009", "011", "013", "0015")
  filter(contains_bcr == TRUE) %>% 
  group_by(cell_type, tissue, sample, clonotype_id_bcr) %>%    # tissue, time
  tally() %>% 
  rename(freq = "n")

# 1) calculate the Gini coefficient per cell_type & per sample
gini_res <- clons_for_gini %>% 
  group_by(cell_type, tissue, sample) %>%        ## tissue, time,
  summarise(across(-clonotype_id_bcr, Gini, na.rm = TRUE)) %>% 
  filter(!freq %in% NaN)

# 1) order the labels
# 2) specify the consistent colors
gini_res$cell_type = factor(gini_res$cell_type, levels=c('B Transitional', 'B Naive',  'B Atypical', 'B Switched', 'Plasmablast'))  #
# gini_res$time = factor(gini_res$time, levels=c('BSE', time_filt))    
gini_res$tissue = factor(gini_res$tissue, levels=c('PBMC', 'CSF'))

# 1) create a column "subject"
gini_res$sample <- as.character(gini_res$sample)
gini_res['subject'] <- sapply(gini_res$sample, function(x){
  strsplit(x, split = "_", fixed = T)[[1]][1]})

gini_res <- gini_res %>% 
  filter(!cell_type %in% c('B Transitional'))
# 1) exclude CD8 Trm because they are only in CSF
# gini_res <- gini_res %>% filter(!cell_type %in% "CD8 Trm")
# gini_res <- gini_res %>% filter(!cell_type %in% "CD8 NKT-like")   # for CSF: W5_10

# 1) run tukey test (NORMAL)
tukey_stat <- gini_res %>%
  group_by(cell_type) %>%
  tukey_hsd(freq ~ tissue, p.adjust.method = 'fdr') %>%     # tissue; time
  filter(group1 == 'BSE') %>%          # PBMC; BSE
  add_xy_position(x = "cell_type") 

# 1) specify your pallete
# BSE: #ea698b, W5_10: #dec9e9, YR1: #a06cd5, YR2: #47126b
my_palette <- c('#ea698b', '#a06cd5')    

# 1) BOXPLOTS WITH STATS: plot the values of Gini coefficient per cell type
fig <- ggplot(gini_res, aes(x = cell_type, y = freq)) +
  geom_boxplot(outlier.shape = NA, aes(fill = tissue), alpha = 0.5) +    # tissue; alpha = 0.5 (0.7 for time)
  geom_point(position=position_dodge(width=0.75), alpha = 0.5,        # alpha = 0.5
             aes(fill = tissue), colour="black", pch=21) +             # tissue
  geom_line(aes(x = as.numeric(factor(cell_type)) + .75 / 4 * ifelse(tissue == "PBMC", -1, 1),   # tissue == "PBMC"; time == "BSE"
                group = interaction(subject, cell_type)), color = "black", linewidth = 0.5, alpha = 0.5) +
  ylab("Gini coefficient") +
  theme_light() +
  theme(axis.title.x = element_blank(),
        axis.title.y=element_text(size=10),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        legend.position="right",
        legend.title=element_text(size=9),
        legend.text=element_text(size=8),
        plot.title = element_text(size = 10, hjust = 0.5))   +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  stat_pvalue_manual(tukey_stat, hide.ns = T) +
  scale_fill_discrete(name = "Tissue") +
  # scale_fill_manual(name = "Time", values = my_palette) +
  ggtitle(paste0("BCR repertoire diversity in BSE: PBMC vs CSF"))   #"BCR repertoire diversity in ", org, ": BSE vs ", time_filt 

# 1) save in this to keep the good asterisk shape
ggsave(filename=paste0('output/figures/BCR/Gini_coef/total/BSE_PBMC_CSF.pdf'), 
       width = 5, height = 3, dpi = 300)
png(filename=paste0('output/figures/BCR/Gini_coef/total/', org, '_BSE_', time_filt, '.png'), 
    width = 5, height = 3, units = "in", res = 300)   # , org, '_BSE_', time_filt, '.png'
plot(fig)
dev.off()

#------------------------------MUTATION FREQ: IMMCANT: CALCULATE #MUTATIONS-----------------

# ## redundant part which descreses the #cells
# # 1) upload *db-pass.tsv and Human_IMGT_IGHV.fasta files
# setwd("/storage1/fs1/martyomov/Active/collaborations/romans/Greg_wu/CLOCK/data/BCR/immcantation/changeo_10x")
# db <- readChangeoDb("1307_0015_BSE_C/filtered_contig_db-pass.tsv")
# ighv <- readIgFasta("IMGT_db/Human_IGHV.fasta")
# 
# # 1) build subject-specific genotype & its germline database
# gt <- inferGenotype(db, germline_db=ighv)
# gtseq <- genotypeFasta(genotype = gt, germline_db = ighv) 
# writeFasta(gtseq, "1307_0015_BSE_C/IGHV_genotype.fasta")


# 1) BASH: run in terminal createGermline.py within each subject-specific file of clonotypes
# 2) merge all merged_p_germ-pass.tsv into one
merged_df <- list.files(pattern = 'merged_p_germ-pass\\.tsv', recursive = TRUE) %>% 
  lapply(function(lst_file) {
    contig_df <- read.csv(lst_file, sep = "\t")
    return(contig_df)
  }) %>%
  bind_rows

path_4 <- '/storage1/fs1/martyomov/Active/collaborations/romans/Greg_wu/CLOCK/output/output/BCR/'
write.csv(merged_df, paste0(path_4, 'immcant_unfilt_clons_subj_germ.csv'))

# 1) load files
path_5 <- '/storage1/fs1/martyomov/Active/collaborations/romans/Greg_wu/CLOCK/'
setwd(path_5)
immc_clons_germ_unf <- read.csv('output/output/BCR/immc_H_subj_clons_germ_unf.csv')
all_contigs_filt <- read.csv('output/output/BCR/all_filtered_BCR.csv')

immc_clons_germ_f <- immc_clons_germ_unf %>% 
  filter(barcode %in% all_contigs_filt$barcode)
write.csv(immc_clons_germ_f, 'output/output/BCR/immc_H_subj_clons_germ_f.csv')

# 1) create a 2nd expData
expData@meta.data <- expData@meta.data %>% 
  select(-c("clonotype_id", "isotype",                  
            "contains_bcr", "sequence_alignment", "germline_alignment_d_mask"))

# 1) leave only "barcode", "sequence_alignment" & "germline_alignment_d_mask" column in immc_filt_clons
# 2) map "sequence_alignment" & "germline_alignment_d_mask" column to GEX
# 3) mark the row order to reverse it after merge function
immc_clons_germ_f_s <- immc_clons_germ_f %>% select(barcode, clone_id, c_call,
                                                    sequence_alignment, germline_alignment_d_mask)
expData@meta.data <- merge(expData@meta.data, immc_clons_germ_f_s, 
                           by=c("barcode"), all.x=TRUE)
rownames(expData@meta.data) <- expData@meta.data$barcode
expData@meta.data <- expData@meta.data[order(expData@meta.data$ordered_id), ]
match.index = match(Cells(expData), immc_clons_germ_f$barcode)
expData@meta.data$contains_bcr <- !is.na(match.index)
expData@meta.data <- expData@meta.data %>% rename(isotype = c_call, clonotype_id = clone_id)

# 1) calculate mutation frequencies
db_mut <- observedMutations(expData@meta.data, 
                            sequenceColumn="sequence_alignment",      # input seq
                            germlineColumn="germline_alignment_d_mask",   # reference seq
                            regionDefinition=NULL,
                            frequency=TRUE,
                            combine=TRUE,
                            nproc=8)
bcr_gex_mut <- db_mut %>%
  filter(contains_bcr %in% TRUE) %>% 
  select(cell_type, time, tissue, mu_freq, isotype,  barcode, sample)
write.csv(bcr_gex_mut, 'output/output/BCR/bcr_mut_freq.csv')
bcr_gex_mut <- read.csv('output/output/BCR/bcr_mut_freq.csv')

# ACROSS CELL TYPES
# 1) define the palette & fix the order of labels 
my_pal <- c('#684756', '#b69121', '#96705b', '#b76935', '#e7bc91')
bcr_gex_mut$cell_type = factor(bcr_gex_mut$cell_type,
                               levels=c('B Transitional', 'B Naive', 'B Atypical', 'B Switched', 'Plasmablast'))

ggplot(bcr_gex_mut, aes(x = cell_type, y = mu_freq, fill = cell_type)) +
  geom_boxplot(outlier.size = 1) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1) +    #, size = 1
  scale_fill_manual(values = my_pal, name = "Cell Type") +
  ggtitle("Mutation frequency by cell type") +
  labs(y = "Mutation frequency") +
  theme_bw() +
  theme(plot.title = element_text(size = 13, hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 10),
        legend.position="none") +
  scale_x_discrete(guide = guide_axis(angle = 45))

ggsave('output/figures/BCR/Mutat_freq/cell_type.png', width = 6, height = 5)
ggsave('output/figures/BCR/Mutat_freq/cell_type.pdf', width = 6, height = 4)

# ACROSS ISOTYPES IN BSE
bcr_gex_mut <- bcr_gex_mut %>% 
  filter(time %in% "BSE")

# 1) define the palette & fix the order of labels 
bcr_gex_mut$isotype = factor(bcr_gex_mut$isotype,
                             levels=c("IGHM", "IGHD", "IGHA1", "IGHA2",
                                      "IGHG1", "IGHG2", "IGHG3", "IGHG4"))
# 1) add custom palette
my_palette <- c('#8ecae6', '#118ab2', '#ffc971', '#ff9505',
                '#ffa69e', '#ff4d6d', '#c9184a', '#590d22')

ggplot(bcr_gex_mut, aes(x = isotype, y = mu_freq, fill = isotype)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1) +    #, size = 1
  scale_fill_manual(values = my_palette, name = "Isotype") +
  ggtitle("Mutation frequency by isotype") +
  labs(y = "Mutation frequency") +
  theme_bw() +
  theme(plot.title = element_text(size = 13, hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 10),
        legend.position="none") +
  scale_x_discrete(guide = guide_axis(angle = 45))

ggsave('output/figures/BCR/Mutat_freq/isotype.png', width = 6, height = 5)
ggsave('output/figures/BCR/Mutat_freq/isotype.pdf', width = 6, height = 4)

#------------------------MUTATION FREQ ACROSS ORGANS/ISOTYPES: ADVANCED BOXPLOT--------------------------------

path_5 <- '/storage1/fs1/martyomov/Active/collaborations/romans/Greg_wu/CLOCK/'
setwd(path_5)
bcr_gex_mut <- read.csv('output/output/BCR/bcr_mut_freq.csv')

mut_res <- bcr_gex_mut %>% 
  filter(time %in% "BSE" & !(isotype %in% c("IGHA2", "IGHG4"))) %>% 
  mutate(time = recode(time, 'W05' = 'W5_10', 'W10' = 'W5_10'))        

# 1) run tukey test (NORMAL)
tukey_stat <- mut_res %>%
  group_by(isotype) %>%
  tukey_hsd(mu_freq ~ tissue, p.adjust.method = 'fdr') %>%     # tissue; time
  filter(group1 == 'PBMC') %>%          # PBMC; BSE
  add_xy_position(x = "isotype")     # cell_type; time 

# 1) order the labels
# 2) specify your pallete
mut_res$isotype = factor(mut_res$isotype,
                         levels=c("IGHM", "IGHD", "IGHA1", "IGHG1", "IGHG2", "IGHG3"))
mut_res$tissue = factor(mut_res$tissue, levels=c('PBMC', 'CSF'))

# 1) BOXPLOTS WITH STATS
fig <-ggplot(mut_res, aes(x = isotype, y = mu_freq)) +    #fig <- 
  geom_boxplot(aes(fill = tissue), alpha = 0.6) +    # tissue; time
  geom_point(position=position_dodge(width=0.75), alpha = 0.6,        # alpha = 0.5
             aes(fill = tissue), colour="black", pch=21) +             # tissue; time
  ylab("Mutation frequency") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y=element_text(size=10),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        legend.position="right",
        legend.title=element_text(size=9),
        legend.text=element_text(size=8),
        plot.title = element_text(size = 10, hjust = 0.5))   +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  stat_pvalue_manual(tukey_stat, hide.ns = T) +
  scale_fill_discrete(name = "Isotype") +
  ggtitle(paste0("Mutation frequency in BSE: PBMC vs CSF")) 

# 1) save in this to keep the good asterisk shape
ggsave(filename=paste0('output/figures/BCR/Mutat_freq/isotype_BSE_PBMC_CSF.pdf'),   
       width = 8, height = 4, dpi = 300)
png(filename=paste0('output/figures/BCR/Mutat_freq/isotype_BSE_PBMC_CSF.png'), 
    width = 8, height = 4, units = "in", res = 300)  
plot(fig)
dev.off()

#------------------------MUTATION FREQ ACROSS ORGANS/CELL TYPES: ADVANCED BOXPLOT--------------------------------

path_5 <- '/storage1/fs1/martyomov/Active/collaborations/romans/Greg_wu/CLOCK/'
setwd(path_5)
bcr_gex_mut <- read.csv('output/output/BCR/bcr_mut_freq.csv')

mut_res <- bcr_gex_mut %>% 
  filter(time %in% "BSE") %>% 
  mutate(time = recode(time, 'W05' = 'W5_10', 'W10' = 'W5_10')) %>%   #'YR1' = 'YR1_2', 'YR2' = 'YR1_2',
  filter(cell_type %in% c('B Transitional', 'B Naive', 'B Atypical', 'B Switched', 'Plasmablast'))        

# 1) run tukey test (NORMAL)
tukey_stat <- mut_res %>%
  group_by(cell_type) %>%
  tukey_hsd(mu_freq ~ tissue, p.adjust.method = 'fdr') %>%     # tissue; time
  filter(group1 == 'PBMC') %>%          # PBMC; BSE
  add_xy_position(x = "cell_type")     # cell_type; time 

# 1) order the labels
# 2) specify your pallete
mut_res$cell_type = factor(mut_res$cell_type, levels=c('B Transitional', 'B Naive', 'B Atypical', 'B Switched', 'Plasmablast'))  
mut_res$tissue = factor(mut_res$tissue, levels=c('PBMC', 'CSF'))
# my_palette <- c('#ea698b', '#dec9e9', '#a06cd5', '#47126b')

# 1) BOXPLOTS WITH STATS
fig <-ggplot(mut_res, aes(x = cell_type, y = mu_freq)) +    #fig <- 
  geom_boxplot(aes(fill = tissue), alpha = 0.6) +    # tissue; time
  geom_point(position=position_dodge(width=0.75), alpha = 0.6,        # alpha = 0.5
             aes(fill = tissue), colour="black", pch=21) +             # tissue; time
  ylab("Mutation frequency") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y=element_text(size=10),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        legend.position="right",
        legend.title=element_text(size=9),
        legend.text=element_text(size=8),
        plot.title = element_text(size = 10, hjust = 0.5))   +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  stat_pvalue_manual(tukey_stat, hide.ns = T) +
  scale_fill_discrete(name = "Tissue") +
  ggtitle(paste0("Mutation frequency in BSE: PBMC vs CSF")) 

# 1) save in this to keep the good asterisk shape
ggsave(filename=paste0('output/figures/BCR/Mutat_freq/cell_type_BSE_PBMC_CSF.pdf'),   
       width = 7, height = 4, dpi = 300)
png(filename=paste0('output/figures/BCR/Mutat_freq/cell_type_BSE_PBMC_CSF.png'), 
    width = 7, height = 4, units = "in", res = 300)  
plot(fig)
dev.off()

#------------------------MUTATION FREQ ACROSS ORGANS/CELL TYPES within IGHM: ADVANCED BOXPLOT--------------------------------

path_5 <- '/storage1/fs1/martyomov/Active/collaborations/romans/Greg_wu/CLOCK/'
setwd(path_5)
bcr_gex_mut <- read.csv('output/output/BCR/bcr_mut_freq.csv')

mut_res <- bcr_gex_mut %>% 
  filter(time %in% "BSE" & isotype %in% "IGHM") %>% 
  mutate(time = recode(time, 'W05' = 'W5_10', 'W10' = 'W5_10')) %>%   #'YR1' = 'YR1_2', 'YR2' = 'YR1_2',
  filter(cell_type %in% c('B Transitional', 'B Naive', 'B Atypical', 'B Switched', 'Plasmablast'))        

# 1) run tukey test (NORMAL)
tukey_stat <- mut_res %>%
  group_by(cell_type) %>%
  tukey_hsd(mu_freq ~ tissue, p.adjust.method = 'fdr') %>%     # tissue; time
  filter(group1 == 'PBMC') %>%          # PBMC; BSE
  add_xy_position(x = "cell_type")     # cell_type; time 

# 1) order the labels
# 2) specify your pallete
mut_res$cell_type = factor(mut_res$cell_type, levels=c('B Transitional', 'B Naive', 'B Atypical', 'B Switched', 'Plasmablast'))  
mut_res$tissue = factor(mut_res$tissue, levels=c('PBMC', 'CSF'))
# my_palette <- c('#ea698b', '#dec9e9', '#a06cd5', '#47126b')

# 1) BOXPLOTS WITH STATS
fig <-ggplot(mut_res, aes(x = cell_type, y = mu_freq)) +    #fig <- 
  geom_boxplot(aes(fill = tissue), alpha = 0.6) +    # tissue; time
  geom_point(position=position_dodge(width=0.75), alpha = 0.6,        # alpha = 0.5
             aes(fill = tissue), colour="black", pch=21) +             # tissue; time
  ylab("Mutation frequency") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y=element_text(size=10),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        legend.position="right",
        legend.title=element_text(size=9),
        legend.text=element_text(size=8),
        plot.title = element_text(size = 10, hjust = 0.5))   +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  stat_pvalue_manual(tukey_stat, hide.ns = T) +
  scale_fill_discrete(name = "Tissue") +
  ggtitle(paste0("Mutation frequency in IGHM of BSE: PBMC vs CSF")) 

# 1) save in this to keep the good asterisk shape
ggsave(filename=paste0('output/figures/BCR/Mutat_freq/IGHM_BSE_PBMC_CSF.pdf'),   
       width = 7, height = 4, dpi = 300)
png(filename=paste0('output/figures/BCR/Mutat_freq/IGHM_BSE_PBMC_CSF.png'), 
    width = 7, height = 4, units = "in", res = 300)  
plot(fig)
dev.off()

#------------------------MUTATION FREQ ACROSS TIME: ADVANCED BOXPLOT--------------------------------

bcr_gex_mut <- read.csv('output/output/BCR/bcr_mut_freq.csv')
org <- 'CSF'

mut_res <- bcr_gex_mut %>% 
  mutate(time = recode(time, 'W05' = 'W5_10', 'W10' = 'W5_10')) %>%   #'YR1' = 'YR1_2', 'YR2' = 'YR1_2',
  filter(tissue %in% org &             
           cell_type %in% c('B Naive', 'B Transitional', 'B Atypical', 'B Switched', 'Plasmablast') &
           time %in% c('BSE', 'W5_10', 'YR1'))        # time %in% c('BSE', 'W5_10', 'YR1') 

# 1) run tukey test (NORMAL)
tukey_stat <- mut_res %>%
  group_by(cell_type) %>%
  tukey_hsd(mu_freq ~ time, p.adjust.method = 'fdr') %>%     # tissue; time
  filter(group1 == 'BSE') %>%          # PBMC; BSE
  add_xy_position(x = "time")     # time

# 1) order the labels
# 2) specify your pallete
mut_res$cell_type <- factor(mut_res$cell_type, levels=c('B Transitional', 'B Naive', 'B Atypical', 'B Switched', 'Plasmablast'))
mut_res$time <- factor(mut_res$time, levels=c('BSE', 'W5_10', 'YR1')) #, 'YR2'
# mut_res$tissue = factor(mut_res$tissue, levels=c('PBMC', 'CSF'))
my_palette <- c('#ea698b', '#dec9e9', '#a06cd5', '#47126b')    # BSE: #ea698b, W5_10: #dec9e9, YR1: #a06cd5, YR2: #47126b

# 1) BOXPLOTS WITH STATS: plot the values of Gini coefficient per cell type
fig <- ggplot(mut_res, aes(x = time, y = mu_freq)) + 
  geom_boxplot(outlier.shape = NA, aes(fill = time), alpha = 0.6) +
  geom_point(position=position_dodge(width=0.75), alpha = 0.6,       
             aes(fill = time), colour="black", pch=21) +             # tissue
  ylab("Mutation frequency") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y=element_text(size=10),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        legend.position="right",
        legend.title=element_text(size=9),
        legend.text=element_text(size=8),
        plot.title = element_text(size = 10, hjust = 0.5))   +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  stat_pvalue_manual(tukey_stat, hide.ns = T) +
  scale_fill_manual(name = "Time", values = my_palette) +
  ggtitle(paste0("Mutation frequency in PBMC across time points groups")) +
  facet_grid(~factor(cell_type, levels=c('B Transitional', 'B Naive', 'B Atypical', 'B Switched', 'Plasmablast')), 
             scales = "free")

# 1) save in this to keep the good asterisk shape
ggsave(filename=paste0('output/figures/BCR/Mutat_freq/PBMC_all_time.pdf'),    #'output/figures/BCR/Mutat_freq/', org, '_BSE_', time_filt, '.pdf'
       width = 10, height = 3.5, dpi = 300)
png(filename=paste0('output/figures/BCR/Mutat_freq/PBMC_all_time.png'), 
    width = 10, height = 3.5, units = "in", res = 300)   # , org, '_BSE_', time_filt, '.png'
plot(fig)
dev.off()

#---------------------------optional: UMAP and mutation freq------------------------------
path <- "/storage1/fs1/martyomov/Active/collaborations/romans/Greg_wu/CLOCK/"
setwd(path)
expData <- readRDS('output/output/scRNA/combined_annot.rds')
bcr_gex_mut <- read.csv('output/output/BCR/bcr_mut_freq.csv')
bcr_gex_mut <- bcr_gex_mut %>% select(barcode, mu_freq)

# 1) add mutation freq column
expData@meta.data <- merge(expData@meta.data, bcr_gex_mut, 
                           by=c("barcode"), all.x=TRUE)
rownames(expData@meta.data) <- expData@meta.data$barcode
expData@meta.data <- expData@meta.data[order(expData@meta.data$ordered_id), ]

# 1) plot mut_freq
FeaturePlot(expData, c('mu_freq')) +
  labs(color = "Substition rate") +
  ggtitle('Mutation frequency in CDR3-H chain region of BCRs') +
  theme(plot.title = element_text(size = 13, hjust = 0.5, face = "plain"),
        legend.title=element_text(size=12))

# 1) plot isotype distribution
expData@meta.data$isotype <- factor(expData@meta.data$isotype,
                                    levels = c("IGHM", "IGHD", "IGHA1", "IGHA2", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "NA"))
isotype_palette = c("IGHM" = "#8ecae6", "IGHD" = "#118ab2", "IGHA1" = "#ffc971", "IGHA2" = "#ff9505",
                    "IGHG1" = "#ffa69e", "IGHG2" = "#ff4d6d", "IGHG3" = "#c9184a", "IGHG4" = "#590d22")
DimPlot(expData, reduction = 'umap', group.by = 'isotype',
        cols = isotype_palette)+ ggtitle(NULL)

ggsave('output/figures/BCR/UMAP/isotype.png', width = 6, height = 4)
ggsave('output/figures/BCR/UMAP/isotype.pdf', width = 6, height = 4)

#---------------------------B CELL LINEAGE TREEs---------------------------

# 1) upload files
path_5 <- '/storage1/fs1/martyomov/Active/collaborations/romans/Greg_wu/CLOCK/'
setwd(path_5)
clons_germ_f <- read.csv('output/output/BCR/immc_H_subj_clons_germ_f.csv')
# expData <- readRDS('output/output/scRNA/combined_annot.rds')
# 
# # 1) add cell_type column to filtered BCR table & save it
# expData_meta_s <- expData@meta.data %>% 
#   filter(contains_bcr %in% TRUE) %>% 
#   select(barcode, cell_type) 
# clons_germ_f <- merge(clons_germ_f, expData_meta_s, by=c("barcode"), all.x=TRUE)
# 
# # 1) modify meta a little bit & save 
# clons_germ_f <- clons_germ_f %>% 
#   mutate(time = recode(time, 'W05' = 'W5_10', 'W10' = 'W5_10')) %>% 
#   rename(isotype = c_call) %>% 
#   mutate(clone_id = paste0('subj_', patient_id, '_clone_', clone_id)) %>% 
#   mutate(cell_type = as.character(cell_type))
# write.csv(clons_germ_f, 'output/output/BCR/immc_H_subj_clons_germ_f.csv')

# 1) look at PBMC samples from subj_6 (BSE, YR1, YR2)
clons_subj_13 <- clons_germ_f %>% 
  filter(patient_id %in% '13')   

# 1) look at clone "478" which wasn't included downstream
freq <- clons_subj_13$clone_id %>% table() %>% as.data.frame()
clons_subj_13_s <- clons_subj_13 %>%
  select('clone_id', 'tissue', 'time', 'cell_type', 'isotype')

# 1) reformat clones with dowser (each row should be a clone)
# it keeps clones with only at least one distinct feature specified in “traits”
clones <- formatClones(clons_subj_13, 
                       traits = c('time', 'cell_type', 'isotype', 'tissue'),   #'time', 'cell_type', 'isotype',
                       minseq = 2, collapse = F,
                       chain = "H", heavy = "IGH", nproc = 8)

# 1) Method 1: pratchet (phangorn): Maximum parsimony-1
trees <- getTrees(clones, build = "pratchet", nproc = 8)
head(trees)

# 1) Method 2: dnapars (PHYLIP): Maximum parsimony-2
dnapars_loc <- '/home/romans/analysis/extra_libs/phylip-3.697/exe/dnapars'
trees = getTrees(clones, build="dnapars", exec=dnapars_loc, nproc = 8)
head(trees)

# 1)  Method 3: pml (phangorn): Standard maximum likelihood
# igphyml_location <- '/home/romans/analysis/extra_libs/igphyml/src/igphyml'
# trees = getTrees(clones, build="pml", trait='tissue', 
#                  igphyml=igphyml_location, nproc = 8)
trees = getTrees(clones, build="pml", nproc = 8)
head(trees)

# 1) Method 4: dnaml (PHYLIP): Standard maximum likelihood
dnaml_loc <- '/home/romans/analysis/extra_libs/phylip-3.697/exe/dnaml'
trees = getTrees(clones, build="dnaml", exec=dnaml_loc, nproc = 8)
head(trees)

# 1) specify a named palette vector
cell_type_palette = c('B Naive' = '#b69121', 'B Transitional' = '#684756', 'B Atypical' = '#96705b',
                      'B Switched' = '#b76935', 'Plasmablast' = '#e7bc91')  #, 'Germline' = '#6c757d'
tissue_palette = c('PBMC' = '#E41A1C', 'CSF' = '#377EB8')
isotype_palette = c("IGHM" = "#8ecae6", "IGHD" = "#118ab2", "IGHA1" = "#ffc971", "IGHA2" = "#ff9505",
                    "IGHG1" = "#ffa69e", "IGHG2" = "#ff4d6d", "IGHG3" = "#c9184a", "IGHG4" = "#590d22")

# 1) to look at the structure of "data" or "tree" obj:
# trees$data[[1]]@data$cell_type
# trees$trees[[1]]$tip.label

# 1) this replacement'll give you NA
clon_id <- 3
trees$trees[[clon_id]]$tip.label <- gsub('Germline', '', trees$trees[[clon_id]]$tip.label)

# 1) Scale branches to mutations rather than mutations/site
# 2) Make fancy tree plot of second largest tree
# trees <- scaleBranches(trees)
plotTrees(trees, scale = 0.03)[[clon_id]] +  #, scale = 1, tips = "cell_type", palette = my_tip_palette
  geom_tippoint(aes(fill = tissue, size = time, color = isotype),
                shape = 21) +   #color = cell_type,
  geom_tiplab(aes(label = cell_type), offset = 0.003) +
  scale_color_manual(values = isotype_palette) +
  scale_fill_manual(values = tissue_palette) +
  theme(plot.title = element_text(size = 15, hjust = 0.5),
        legend.position="right",
        legend.title=element_text(size=12),
        legend.text=element_text(size=10)) +
  labs(fill="Tissue", size="Time", color="Isotype") +
  guides(fill = guide_legend(order = 1),
         size = guide_legend(order = 2),
         color = guide_legend(order = 3)) +
  ggtitle(paste0('B cell tree of ', clones$clone_id[[clon_id]])) 

ggsave(paste0('output/figures/BCR/Phyl_trees/', clones$clone_id[[clon_id]], '_new.pdf'),
       width = 7, height = 6.5)
ggsave(paste0('output/figures/BCR/Phyl_trees/', clones$clone_id[[clon_id]], '_new.png'), 
       width = 7, height = 6.5)

#-----------------------------Follow-ups on Brian's questions--------------------------

path <- "/storage1/fs1/martyomov/Active/collaborations/romans/Greg_wu/CLOCK/"
setwd(path)
expData <- readRDS('output/output/scRNA/combined_annot.rds')
bcr_gex_mut <- read.csv('output/output/BCR/bcr_mut_freq.csv')
bcr_gex_mut <- bcr_gex_mut %>% select(barcode, mu_freq)

# 1) add mutation freq column
expData@meta.data <- merge(expData@meta.data, bcr_gex_mut, 
                           by=c("barcode"), all.x=TRUE)
rownames(expData@meta.data) <- expData@meta.data$barcode
expData@meta.data <- expData@meta.data[order(expData@meta.data$ordered_id), ]

# 1) subset the expData
b_naive_CSF <- subset(expData, (cell_type %in% 'B Naive' & tissue %in% 'CSF' & contains_bcr %in% TRUE))

DimPlot(b_naive_CSF, reduction = 'umap', group.by = 'cell_type')
FeaturePlot(expData, c('mu_freq')) +
  labs(color = "Substition rate") +
  ggtitle('Mutation frequency in CDR3-H chain region of BCRs') +
  theme(plot.title = element_text(size = 13, hjust = 0.5, face = "plain"),
        legend.title=element_text(size=12))

ggsave('output/figures/BCR/UMAP/mut_freq.png', width = 6, height = 4)
ggsave('output/figures/BCR/UMAP/mut_freq.pdf', width = 6, height = 4)













