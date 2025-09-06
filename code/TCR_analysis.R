library(dplyr)
library(ggplot2)
library(stringr)
library(Seurat)
library(forcats)
library(mgsub)
library(DescTools)
library(ggpubr)
library(rstatix)
library(tidyverse)
library(ggbreak) 
library(scRepertoire)
setwd("/storage1/fs1/martyomov/Active/collaborations/romans/Greg_wu/CLOCK")
#

# -----------------------------Combine all TCRs into one dataframe--------------

path_2 <- "/storage1/fs1/martyomov/Active/collaborations/romans/Greg_wu/CLOCK/data/TCR"
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

path_3 <- '/storage1/fs1/martyomov/Active/collaborations/romans/Greg_wu/CLOCK/output/output/TCR/'
write.csv(merged_df, paste0(path_3, 'all_unfiltered_TCR.csv'))

#----------------------------Filter out cell with "non-canonical" TCRs----------

# 1) set env & load the data
path <- "/storage1/fs1/martyomov/Active/collaborations/romans/Greg_wu/CLOCK/"
setwd(path)
all_contigs <- read.csv(paste0(path, "output/output/TCR/all_unfiltered_TCR.csv"))

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


# 1) count #chains & leave only barcodes with 1 alpha: 1 betta
all_contigs_filt$n_alpha <- sapply(1:nrow(all_contigs_filt), function(i) str_count(all_contigs_filt$cdr3s_aa[i], "TRA:"))
all_contigs_filt$n_beta <- sapply(1:nrow(all_contigs_filt), function(i) str_count(all_contigs_filt$cdr3s_aa[i], "TRB:"))
all_contigs_filt <- all_contigs_filt %>% filter(n_alpha == 1 & n_beta == 1)

# 1) extract cell IDs from GEX
# 2) remove TCRs that are not in GEX
expData <- readRDS('output/output/scRNA/combined_annot.rds')
tcr_present <- rownames(expData@meta.data)
all_contigs_filt <- all_contigs_filt %>% 
  filter(barcode %in% tcr_present)

# 1) save the checkpoint & load it again
write.csv(all_contigs_filt, 'output/output/TCR/all_filtered_TCR.csv')
all_contigs_filt <- read.csv('output/output/TCR/all_filtered_TCR.csv')
expData <- readRDS('output/output/scRNA/combined_annot.rds')

# 1) add all columns from TCR data to Seurat obj
# 2) mark the row order to reverse it after merge function
expData@meta.data <- expData@meta.data %>% 
  select(-barcode) %>% rename(barcode = unique_barcode)
expData@meta.data$ordered_id  <- 1:nrow(expData@meta.data)
expData@meta.data <- merge(expData@meta.data, all_contigs_filt, 
                           by=c("barcode"), all.x=TRUE)
rownames(expData@meta.data) <- expData@meta.data$barcode
expData@meta.data <- expData@meta.data[order(expData@meta.data$ordered_id), ]

# 1) add TCR data to Seurat obj
# 2) filter out non-T cells which have TCRs & save GEX
match.index = match(Cells(expData), all_contigs_filt$barcode)
expData@meta.data$contains_tcr <- !is.na(match.index)
# expData@meta.data <- expData@meta.data %>%
#   filter(!(cell_type %in% c('NK CD56dim','NK CD56bright', 'ILC', 'CD14 Mono', 'CD16 Mono', 'cDC', 'pDC', 'Microglia-like', 'B Naive', 'B Transitional', 'B Atypical', 'B Switched', 'Plasmablast', 'gdT2')
#            & contains_tcr == TRUE))
# cellsToKeep <- rownames(expData@meta.data)
# write.csv(cellsToKeep, 'output/output/scRNA/cellsToKeep.csv')
# expData <- expData[, colnames(expData) %in% cellsToKeep]
# saveRDS(expData, 'output/output/scRNA/combined_annot.rds')

# 1) visualize cells with TCRs on the UMAP
highlighted_cells <- Cells(expData)[which(expData$contains_tcr)]
Idents(expData) <- expData$cell_type
DimPlot(expData, reduction = 'umap', cells.highlight = highlighted_cells, label.size = 2.5,
        repel=TRUE, label = TRUE, cols = 'grey', pt.size = 0.1, sizes.highlight = 0.5) +
  scale_color_manual(labels = c("Non-clonal", "Clonal"),
                     values = c("grey", "#d90429")) 
ggsave('output/ <- ures/TCR/UMAP/tcr_umap.png', width = 10, height = 7)
ggsave('output/figures/TCR/UMAP/tcr_umap.pdf', width = 7, height = 5)

# 1) recalculate TCR clonotypes based on cdr3s_aa
all_contigs_filt <- all_contigs_filt %>% 
  group_by(cdr3s_aa) %>% 
  mutate(clonotype_id = cur_group_id()) %>% 
  mutate(clonotype_id = paste0('clonotype_', clonotype_id))

# 1) add high level annotation column for CD4 & CD8 T cells
expData@meta.data <- expData@meta.data %>% 
  mutate(cell_type_2 = gsub("CD8 Naive|CD8 Tem GZMK|CD8 Tcm CCR4|CD8 Tem GZMB|Proliferative T/NK|CD8 Tcm CCR4+|CD8 NKT-like|CD8 Trm|CD8 Temra", "CD8 T cell", 
                            gsub("CD4 Th1|CD4 Th17|CD4 Naive|CD4 Th22|CD4 Treg naive|CD4 Treg memory|CD4 Temra|CD4 exhausted|Proliferative T/NK", "CD4 T cell",
                                 cell_type))) %>% 
  mutate(cell_type_2 = recode(cell_type_2, 'CD8 T cell+' = 'CD8 T cell',
                              'CD8 T cell-' = 'CD8 T cell'))


#-------------Dotplot with clonotypes frequencies across all T cells-------------

# 1) count # unique clones for cdr3_aa
contig_short <- all_contigs_filt
uniq_clones <- contig_short$clonotype_id %>% 
  table() %>% as.data.frame()

# 1) compress it, by calculating #clones within a clonotype
uniq_clones_table_2 <- uniq_clones$Freq %>%
  table() %>% as.data.frame()

# 1) rename columns & calculate log(1 + x) to rescale freq
uniq_clones_table_2 <- uniq_clones_table_2 %>% 
  rename(clone_size = '.', freq = Freq) %>% 
  mutate(log1p_frequency = log1p(freq)) 

# 1) plot
ggplot(uniq_clones_table_2, aes(x = fct_inseq(clone_size),y = log1p_frequency)) +
  geom_point() +
  xlab("Clone size") +
  ylab("log1p of #clonotypes") +
  theme(axis.text.x = element_text(angle=45, hjust=1)) 

ggsave('/storage1/fs1/martyomov/Active/collaborations/romans/Greg_wu/CLOCK/output/figures/TCR/Freq_dotplot/all_clones_doplot.png',
       dpi = 300, width = 12, height = 7)
ggsave('/storage1/fs1/martyomov/Active/collaborations/romans/Greg_wu/CLOCK/output/figures/TCR/Freq_dotplot/all_clones_doplot.pdf',
       dpi = 300, width = 12, height = 7)

#-------------Dotplot with subjects' frequencies across clonotypes-------------

# 1) count # unique clones for cdr3_aa
contig_short <- all_contigs_filt
contig_short['subject'] <- sapply(contig_short$barcode, function(x){strsplit(x, split = "_", fixed = T)[[1]][1]})
uniq_clones <- contig_short$clonotype_id %>% 
  table() %>% as.data.frame()

# 1) compress it, by calculating #subjects shared the same clonotype
n_subject <- contig_short %>% 
  group_by(clonotype_id) %>% 
  summarise(n_subject = n_distinct(subject))

# 1) rename columns & merge two df
uniq_clones <- uniq_clones %>% 
  rename(clonotype_id = '.')
clons_subjs <- merge(uniq_clones, n_subject, by = 'clonotype_id')

# 1) plot
ggplot(clons_subjs, aes(x = fct_inseq(as.factor(Freq)), y = n_subject)) +
  geom_point() +
  xlab("Clone size") +
  ylab("# subjects") +
  theme(axis.text.x = element_text(angle=45, hjust=1)) 

ggsave('/storage1/fs1/martyomov/Active/collaborations/romans/Greg_wu/CLOCK/output/figures/TCR/Freq_dotplot/n_subjects_clonotype.png',
       dpi = 300, width = 12, height = 7)
ggsave('/storage1/fs1/martyomov/Active/collaborations/romans/Greg_wu/CLOCK/output/figures/TCR/Freq_dotplot/n_subjects_clonotype.pdf',
       dpi = 300, width = 12, height = 7)

#-------------------------------------DYNAMICS OF INDIVIDUAL CLONOTYPES-------------------

# 1) set env & load the data
path <- "/storage1/fs1/martyomov/Active/collaborations/romans/Greg_wu/CLOCK/"
setwd(path)
expData <- readRDS('output/output/scRNA/combined_annot.rds')
all_contigs <- read.csv("output/output/TCR/filtered_TCR.csv")

# 1) count # unique clones for cdr3_aa
contig_short <- all_contigs
contig_short['subject'] <- sapply(contig_short$barcode, function(x){strsplit(x, split = "_", fixed = T)[[1]][1]})
uniq_clones <- contig_short$clonotype_id %>% 
  table() %>% as.data.frame() %>% 
  rename(clonotype_id = '.', clone_size = Freq)

# 1) subset barcodes of the biggest clonotype (e.g., clonotype_24487)
clonotype_comp <- expData@meta.data %>% 
  filter(clonotype_id_tcr %in% 'clonotype_7449') %>% 
  select(clonotype_id_tcr, barcode, subject, tissue, time, cell_type, cell_type_2)

clonotype_comp <- clonotype_comp %>% filter(!cell_type_2 %in% 'CD4 T cell')

# 1) look at CD4 T cells from one clone on the UMAP
cd4_clone <- clonotype_comp %>% filter(cell_type_2 %in% 'CD4 T cell')
cd4_clone_gex <- subset(expData, (barcode %in% cd4_clone$barcode))
DimPlot(cd4_clone_gex, reduction = 'umap', group.by = 'cell_type')
FeaturePlot(cd4_clone_gex, c('CXCR3', 'GZMH'), ncol = 1)
ggsave('output/figures/TCR/UMAP/mark2_cd4_clone_7449.png', width = 4, height = 7)
ggsave('output/figures/TCR/UMAP/mark2_cd4_clone_7449.pdf', width = 4, height = 7)

# 1) prepare automatic info for the plot title
clonotype_id <- clonotype_comp$clonotype_id_tcr[1]
clon_size <- clonotype_comp %>% nrow()
subj <- clonotype_comp$subject[1]

# 1) count #clones per tissue, time & cell_type
clonotype_comp_gr <- clonotype_comp %>% 
  group_by(tissue, time, cell_type) %>% 
  tally() %>% 
  rename(clone_size = n)

# 1) fix the order of labels
clonotype_comp_gr$tissue <- factor(clonotype_comp_gr$tissue, levels=c("PBMC", "CSF"))
clonotype_comp_gr$time <- factor(clonotype_comp_gr$time, levels=c("BSE", "W05", "W10", "YR1", "YR2"))
clonotype_comp_gr$cell_type <- factor(clonotype_comp_gr$cell_type,
                                      levels=c('CD4 Naive', 'CD4 Treg naive', 'CD4 Treg memory', 'CD4 Th1',
                                               'CD4 Th17', 'CD4 Th22', 'CD4 Temra', 'CD4 exhausted',
                                               'CD8 Naive', 'CD8 Tcm CCR4-', 'CD8 Tcm CCR4+', 'CD8 Trm',
                                               'CD8 Tem GZMK+', 'CD8 Tem GZMB+', 'CD8 Temra', 'CD8 NKT-like',
                                               'gdT2', 'MAIT', 'Proliferative T/NK'))
my_palette <- c('#d60000',  '#ff7266', '#edb8b8', '#e36414',
                '#ff9f1c', '#ffbf69', '#f9cb9cff', '#fdf490', 
                '#a3b18a', '#588157', '#3a5a40', '#bce784',
                '#55a630', '#007f5f', '#a5be00', '#74d3ae',
                '#0774d8', '#00acc6','#9651c4')
names(my_palette) <- levels(clonotype_comp_gr$cell_type)

# 1) plot using stacked barplot
ggplot(clonotype_comp_gr) +
  geom_bar(aes(x = time, y = clone_size, fill = cell_type),
            stat = "identity", color='black') +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 10)) +
  ggtitle(paste0(clonotype_id, " with clone size = ", clon_size, " in subject ", subj)) +
  xlab("Time") +
  ylab("# clones") +
  scale_fill_manual(name = "Cell type", values = my_palette) +
  facet_grid(~tissue)

ggsave(paste0('output/figures/TCR/Freq_barplot/subject_big_clon/subject_', subj, '_clon_size_', clon_size, '.png'),
       width = 6, height = 3.5,  dpi = 300)
ggsave(paste0('output/figures/TCR/Freq_barplot/subject_big_clon/subject_', subj, '_clon_size_', clon_size, '.pdf'),
       width = 6, height = 3.5,  dpi = 300)

#-----------------------SCREP: CLONAL DYNAMICS-----------------------------

# 1) load the data
path <- "/storage1/fs1/martyomov/Active/collaborations/romans/Greg_wu/CLOCK/"
setwd(path)
all_contigs_filt <- read.csv('output/output/TCR/all_filtered_TCR.csv')
all_contigs_unf <- read.csv('output/output/TCR/all_unfiltered_TCR.csv')
expData <- readRDS('output/output/scRNA/combined_annot.rds')

tcr_true <- expData@meta.data %>% filter(contains_tcr %in% TRUE)
new_contigs_filt <- all_contigs_filt %>% 
  filter(barcode %in% tcr_true$barcode)
write.csv(new_contigs_filt, 'output/output/TCR/filtered_TCR.csv')

# 1) modify & filter cells by barcode
contigs_full_f <- all_contigs_unf %>% 
  mutate(barcode = gsub("-1", "", unique_barcode)) %>% 
  filter(barcode %in% new_contigs_filt$barcode)

# 1) reformat data into scRep format of lists
path_2 <- "/storage1/fs1/martyomov/Active/collaborations/romans/Greg_wu/CLOCK/data/TCR"
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
      filter(unique_barcode %in% new_contigs_filt$barcode) %>% 
      select(-unique_barcode)
    return(contig_df)
  })

path_3 <- '/storage1/fs1/martyomov/Active/collaborations/romans/Greg_wu/CLOCK/output/output/TCR/'
saveRDS(merged_df, paste0(path_3, "full_filtered_TCR_list_2.RData"))
full_filtered_TCR <- readRDS(paste0(path_3, 'full_filtered_TCR_list_2.RData'))


for (i in 1:length(full_filtered_TCR)){
  print(unique(full_filtered_TCR[[i]]$sample)) 
}

# 1) combine contigs into clones
combined.TCR <- combineTCR(full_filtered_TCR, 
                           samples = c("001_YR2_CSF" ,"001_YR2_PBMC" ,"0015_BSE_CSF" ,"0015_BSE_PBMC"
                                       ,"002_YR2_PBMC" ,"003_YR2_PBMC" ,"004_YR2_PBMC" ,"005_BSE_CSF"
                                       ,"005_BSE_PBMC" ,"005_YR1_CSF" ,"005_YR1_PBMC" ,"006_BSE_CSF"
                                       ,"006_BSE_PBMC" ,"006_YR1_CSF" ,"006_YR1_PBMC" ,"006_YR2_PBMC"
                                       ,"007_W10_PBMC" ,"007_YR2_PBMC" ,"009_BSE_CSF" ,"009_BSE_PBMC"
                                       ,"009_YR2_PBMC" ,"010_BSE_PBMC" ,"010_YR2_PBMC" ,"011_BSE_CSF"
                                       ,"011_BSE_PBMC" ,"011_W10_CSF" ,"011_W10_PBMC" ,"011_YR2_PBMC"
                                       ,"012_BSE_PBMC" ,"012_YR2_PBMC" ,"013_BSE_CSF" ,"013_BSE_PBMC"
                                       ,"013_W05_CSF" ,"013_W05_PBMC"),
                           removeNA = TRUE, 
                           removeMulti = TRUE, 
                           filterMulti = FALSE)

# 1) manually find top 10 clones at BSE within each sample of interest
bse_df <- as.data.frame(combined.TCR$'005_BSE_CSF')
clonotypes <- bse_df %>% 
  group_by(CTaa) %>% 
  tally() %>% 
  slice_max(order_by = n, n = 10, with_ties = FALSE)          
top_10 <- clonotypes$CTaa

# 1) plot alluvial stacked barplot
# "005_BSE_PBMC",  "005_YR1_PBMC", "005_BSE_CSF", "005_YR1_CSF"
# "006_BSE_PBMC", "006_YR1_PBMC", "006_YR2_PBMC", "006_BSE_CSF", "006_YR1_CSF"
# "011_BSE_PBMC" ,"011_W10_PBMC" ,"011_YR2_PBMC", "011_BSE_CSF" ,"011_W10_CSF"
# "013_BSE_PBMC","013_W05_PBMC", "013_BSE_CSF", "013_W05_CSF" 
p <- clonalCompare(combined.TCR, top.clones = 10, 
                   clones = top_10,
                   samples = c("005_BSE_CSF", "005_YR1_CSF" ), 
                   cloneCall="aa", graph = "alluvial", palette = "inferno") 

p + theme(axis.text.x = element_text(angle=45, hjust=1)) +
  ylab("Proportion of clonotypes") +
  ggtitle("Top 10 clonotypes in CSF of subj_5")

ggsave(paste0(path, 'output/figures/TCR/Clone_barplot/subj_5_CSF.png'), 
       width = 5, height = 4)
ggsave(paste0(path, 'output/figures/TCR/Clone_barplot/subj_5_CSF.pdf'), 
       width = 5, height = 4)

#-------------------------------------UMAP WITH CLONE SIZE-----------------------

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

#-------------------------------------UMAP-MAX CLONE PLOT-----------------------

# 1) UMAP-MAX CLONE: calculate frequncies of clonotypes grouped by tissue & cell_type_2 & pick the biggest
freq_tcr <- expData@meta.data %>%
  filter(contains_tcr == TRUE) %>% 
  group_by(tissue, cell_type_2, clonotype_id_tcr) %>% 
  tally() %>% 
  slice(which.max(n))

# 1) visualize the  on the UMAP
tissue <- 'PBMC'
highlight_clone <- Cells(expData)[which(expData$clonotype_id_tcr == 'clonotype_27369' &
                                          expData$tissue == tissue & 
                                          expData$cell_type_2 == 'CD4 T cell')]
highlight_cd4 <- Cells(expData)[which(expData$tissue == tissue & 
                                        expData$cell_type_2 == 'CD4 T cell' &
                                        !(expData$barcode %in% highlight_clone))]
Idents(expData) <- expData$cell_type
DimPlot(expData, reduction = 'umap', cells.highlight = list(highlight_clone, highlight_cd4),
        repel=TRUE, label = FALSE, cols = 'grey', pt.size = 0.1, sizes.highlight = c(0.7, 0.2),   
        cols.highlight = c(alpha("#8EA2B4", 1), "#d90429")) +
  scale_color_manual(labels = c("Non-clonal", "CD4 T cells", "MS Clone"),   #max CD4 T cell clone in PBMC
                     values = c("grey", "#ff758f", "#d90429")) +        #"#8EA2B4"  OR "#ff758f" -pink
  ggtitle(paste0("Maximum MS Clone 1 in ", tissue)) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(paste0('output/figures/TCR/UMAP/cd4_', tissue, '_max_clone.png'), width = 7, height = 5)  #10*7
ggsave(paste0('output/figures/TCR/UMAP/cd4_', tissue, '_max_clone.pdf'), width = 7, height = 5)

#-----------------------------------------DONUT CHART---------------------------

# 1) DONUT CHART: calculate frequncies of clonotypes grouped by tissue & cell_type_2 & pick the biggest
freq_tcr <- expData@meta.data %>%
  filter(subject %in% c('005', '006', '009', '011', '013', '0015')) %>%   #filter matched samples
  mutate(time_comb = time) %>% 
  mutate(time_comb = recode(time_comb, 'W05' = 'W5_10', 'W10' = 'W5_10')) %>%   #'YR1' = 'YR1_2', 'YR2' = 'YR1_2',
  filter(contains_tcr == TRUE) %>% 
  group_by(cell_type_2, tissue, time_comb, clonotype_id) %>% 
  tally() %>% 
  filter(cell_type_2 %in% c('CD8 T cell', 'CD4 T cell'))     

# 1) keep only those groups in which you are interested in
# organ <- 'PBMC'
# time <- 'YR2'
freq_tcr_filt <- freq_tcr %>% 
  filter(cell_type_2 %in% 'CD8 T cell' &
           tissue %in% c("PBMC", "CSF") &               # organ
           time_comb %in% c('BSE')) %>%   #, time
  rename(clon_size = "n")
freq_tcr_filt$clon_size %>% table()

# 1) count the total #clones per time group
# 2) calculate in % the contribution of clones in each group
num_tot_cells <- freq_tcr_filt %>% 
  group_by(tissue) %>%   #tissue, time_comb
  summarise(total_clon_size = sum(clon_size)) 
freq_tcr_filt$total_clon_size <- with(freq_tcr_filt, 
                                      num_tot_cells$total_clon_size[match(tissue, num_tot_cells$tissue)])   #time_comb, tissue
freq_tcr_filt <- freq_tcr_filt %>% 
  mutate(percent = clon_size/total_clon_size * 100)

# 1) create the colomn "clone_group" using when_case & split into 4 groups (unique, 2-4, 5-9, >=10)
freq_tcr_filt <- freq_tcr_filt %>% 
  mutate(clone_group = case_when(clon_size == 1 ~ 'Unique',        ## for Greg: clonotype_id == 'clonotype_24487' ~ 'clonotype_24487',
                                 clon_size >=2 & clon_size <= 4 ~ '2-4',
                                 clon_size >=5 & clon_size <= 9 ~ '5-9',
                                 clon_size >= 10 ~ '>=10'))

# 1) collapse based on clone_group/Unique group
freq_tcr_filt <- freq_tcr_filt %>% 
  group_by(cell_type_2, tissue, time_comb, 
           grp = case_when(clone_group %in% "Unique" ~ 0, TRUE ~ row_number())) %>% 
  summarise(across(starts_with(c("clon_size", "percent")), sum),
            clone_group = str_c(clone_group, collapse=",")) %>% 
  select(-grp) %>% 
  mutate(clone_group = gsub("^Uni.*", "Unique", clone_group))

# 1) collapse based on "2-4" group
freq_tcr_filt <- freq_tcr_filt %>% 
  group_by(cell_type_2, tissue, time_comb, 
           grp = case_when(clone_group %in% "2-4" ~ 0, TRUE ~ row_number())) %>% 
  summarise(across(starts_with(c("clon_size", "percent")), sum),
            clone_group = str_c(clone_group, collapse=",")) %>% 
  select(-grp) %>% 
  mutate(clone_group = gsub("^2.*", "2-4", clone_group))

# 1) order the labels
freq_tcr_filt$clone_group %>% unique()
freq_tcr_filt$clone_group = factor(freq_tcr_filt$clone_group,
                                   levels=c("Unique", "2-4", "5-9", ">=10"))   #, "clonotype_24487"
freq_tcr_filt$tissue = factor(freq_tcr_filt$tissue, levels=c("PBMC", "CSF"))

# 1) plot the Donut chart
colors <- c("#adb5bd", "#ff8fa3", "#ff4d6d", "#a4133c")   #, "#e09f3e"
ggplot(data = freq_tcr_filt, aes(x=2, y = percent))+    #, fill = clone_group
  geom_bar(stat="identity", color='black')+
  geom_bar(stat="identity", aes(fill = clone_group))+
  # geom_col(color='black')+
  coord_polar("y", start = 0) + 
  facet_wrap(~ tissue) +    #tissue, time_comb
  theme(panel.background = element_blank(),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(), 
        strip.text = element_text(size = 15)) +
  scale_fill_manual(values = colors, name = 'No. of clones') +
  xlim(0.5, 2.5)  
ggsave('output/figures/TCR/Donut_chart/across_matched/CD8_BSE_PBMC_CSF.png', width = 10, height = 7)
ggsave('output/figures/TCR/Donut_chart/across_matched/CD8_BSE_PBMC_CSF.pdf', width = 10, height = 7)
# ggsave(paste0('output/figures/TCR/Donut_chart/across_matched/CD4_', time, '_', organ, '.png'), width = 10, height = 7)
# ggsave(paste0('output/figures/TCR/Donut_chart/across_matched/CD4_', time, '_', organ, '.pdf'), width = 10, height = 7)

#--------------------------------------GINI COEFFICIENT-------------------------

# 1) GINI COEFFICIENT: toy example & intuition behind Gini coef.
x1=c(2,2,2,2)
x2=c(0,1,2,5)
Gini(x1)
Gini(x2)

# load 
expData <- readRDS('output/output/scRNA/combined_annot.rds')

### 1) remove useless columns FOR SIMPLE BOXPLOT
# 2) and count clone sizes
org <- "CSF"
clons_for_gini <- expData@meta.data %>% 
  filter(tissue %in% org &
           cell_type %in% c('CD4 Naive', 'CD4 Treg naive', 'CD4 Treg memory', 'CD4 Th1', 'CD4 Th17', 'CD4 Th22', 'CD4 Temra',
                            'CD4 exhausted', 'CD8 Naive', 'CD8 Tcm CCR4-', 'CD8 Tem GZMK+', 'CD8 Tem GZMB+') ) %>%   
  filter(contains_tcr == TRUE) %>% 
  group_by(cell_type, sample, clonotype_id_tcr) %>%    
  tally() %>% 
  rename(freq = "n")

# 1) calculate the Gini coefficient per cell_type & per sample
gini_res <- clons_for_gini %>% 
  group_by(cell_type, sample) %>%       
  summarise(across(-clonotype_id_tcr, Gini, na.rm = TRUE)) %>% 
  filter(!freq %in% NaN)

# 1) order the labels
# 2) specify the consistent colors
gini_res$cell_type = factor(gini_res$cell_type, levels=c('CD4 Naive', 'CD4 Treg naive', 'CD4 Treg memory', 'CD4 Th1', 'CD4 Th17', 'CD4 Th22', 'CD4 Temra',
                                                         'CD4 exhausted', 'CD8 Naive', 'CD8 Tcm CCR4-', 'CD8 Tem GZMK+', 'CD8 Tem GZMB+'))
my_palette <- c('#d60000',  '#ff7266', '#edb8b8', '#e36414', '#ff9f1c', '#ffbf69', '#f9cb9cff', 
                '#fdf490', '#018700', '#004b00',  '#5d7e66', '#95b577')

## 1) BOXPLOTS WITHOUT STATS: plot the values of Gini coefficient per cell type
ggboxplot(gini_res, x = 'cell_type', y = 'freq', fill = 'cell_type') +
  geom_point(position=position_dodge(width=0.8),
             aes(fill = cell_type), colour="black", pch=21) +
  ylab(paste0("Gini coefficient")) +
  theme_light() +
  theme(axis.title.x = element_blank(),
        axis.title.y=element_text(size=10),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        legend.position="none",
        plot.title = element_text(size = 10, hjust = 0.5))   +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  scale_fill_manual(values = my_palette) +
  ggtitle(paste0("TCR repertoire diversity in ", org))

ggsave(paste0('output/figures/TCR/Gini_coef/CD4_CD8_', org, '_all.png'),
       width = 7, height = 4,  dpi = 300)

### 1) remove useless columns FOR COMPLEX/ADVANCED BOXPLOT
# 2) and count clone sizes
org <- 'PBMC'
time_filt <- 'W5_10'
clons_for_gini <- expData@meta.data %>% 
  mutate(time = recode(time, 'W05' = 'W5_10', 'W10' = 'W5_10')) %>%   #'YR1' = 'YR1_2', 'YR2' = 'YR1_2',
  filter(tissue %in% org &             ##
           cell_type %in% c('CD4 Naive', 'CD4 Treg naive', 'CD4 Treg memory', 'CD4 Th1', 'CD4 Th17', 'CD4 Th22', 'CD4 Temra', 'CD4 exhausted',
                            'CD8 Naive', 'CD8 Tcm CCR4-', 'CD8 Tcm CCR4+', 'CD8 Trm', 'CD8 Tem GZMK+', 'CD8 Tem GZMB+', 'CD8 Temra', 'CD8 NKT-like') & 
           time %in% c("BSE", time_filt)                        ## , time_filt
  ) %>%   ## subject %in% c("005", "006", "009", "011", "013", "0015")
  filter(contains_tcr == TRUE) %>% 
  group_by(cell_type, time, sample, clonotype_id_tcr) %>%    # tissue, time, clonotype_id
  tally() %>% 
  rename(freq = "n")

# 1) calculate the Gini coefficient per cell_type & per sample
gini_res <- clons_for_gini %>% 
  group_by(cell_type, time, sample) %>%        ## tissue, time,
  summarise(across(-clonotype_id_tcr, Gini, na.rm = TRUE)) %>%  #clonotype_id
  filter(!freq %in% NaN)

# 1) order the labels
# 2) specify the consistent colors
gini_res$cell_type = factor(gini_res$cell_type, levels=c('CD4 Naive', 'CD4 Treg naive', 'CD4 Treg memory', 'CD4 Th1', 'CD4 Th17', 'CD4 Th22', 'CD4 Temra', 'CD4 exhausted',
                                                         'CD8 Naive', 'CD8 Tcm CCR4-', 'CD8 Tcm CCR4+', 'CD8 Trm', 'CD8 Tem GZMK+', 'CD8 Tem GZMB+', 'CD8 Temra', 'CD8 NKT-like'))
gini_res$time = factor(gini_res$time, levels=c('BSE', time_filt))    
# gini_res$tissue = factor(gini_res$tissue, levels=c('PBMC', 'CSF'))

# 1) create a column "subject"
gini_res$sample <- as.character(gini_res$sample)
gini_res['subject'] <- sapply(gini_res$sample, function(x){
  strsplit(x, split = "_", fixed = T)[[1]][1]})

# 1) exclude CD8 Trm because they are only in CSF
gini_res <- gini_res %>% filter(!cell_type %in% "CD8 Trm")
# gini_res <- gini_res %>% filter(!cell_type %in% "CD8 NKT-like")   # for CSF: W5_10

# 1) run tukey test (NORMAL)
tukey_stat <- gini_res %>%
  group_by(cell_type) %>%
  tukey_hsd(freq ~ time, p.adjust.method = 'fdr') %>%     # tissue; time
  filter(group1 == 'BSE') %>%          # PBMC; BSE
  add_xy_position(x = "cell_type") 

# 1) specify your pallete
# BSE: #ea698b, W5_10: #dec9e9, YR1: #a06cd5, YR2: #47126b
my_palette <- c('#ea698b', '#dec9e9')    

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
  ggtitle(paste0("TCR repertoire diversity in ", org, ": BSE vs ", time_filt))   #"TCR repertoire diversity in BSE: PBMC vs CSF"

# 1) save in this to keep the good asterisk shape
ggsave(filename=paste0('output/figures/TCR/Gini_coef/total/', org, '_BSE_', time_filt, '.pdf'), 
       width = 10, height = 4, dpi = 300)
png(filename=paste0('output/figures/TCR/Gini_coef/total/', org, '_BSE_', time_filt, '.png'), 
    width = 10, height = 4, units = "in", res = 300)   # , org, '_BSE_', time_filt, '.png'
plot(fig)
dev.off()

### 1) CD4/CD8 T cells: BSE vs W5_10 vs YR1 vs YR2 (for CSF & PBMC) 
# 1) set env & load the data
path <- "/storage1/fs1/martyomov/Active/collaborations/romans/Greg_wu/CLOCK/"
setwd(path)
expData <- readRDS('output/output/scRNA/combined_annot.rds')
all_contigs_filt <- read.csv('output/output/TCR/all_filtered_TCR.csv')

# 1) add "contains_tcr" colm
match.index = match(Cells(expData), all_contigs_filt$barcode)
expData@meta.data$contains_tcr <- !is.na(match.index)

# 1) leave only CD4 & CD8 T cells
# 'CD4 Naive', 'CD4 Treg naive', 'CD4 Treg memory', 'CD4 Th1', 'CD4 Th17', 'CD4 Th22', 'CD4 Temra', 'CD4 exhausted'
clons_for_gini <- expData@meta.data %>% 
  mutate(time = recode(time, 'W05' = 'W5_10', 'W10' = 'W5_10')) %>%  
  filter(cell_type %in% c('CD8 Naive', 'CD8 Tcm CCR4-', 'CD8 Tcm CCR4+', 'CD8 Trm', 'CD8 Tem GZMK+', 'CD8 Tem GZMB+', 'CD8 Temra', 'CD8 NKT-like')) %>%   #
  filter(contains_tcr == TRUE) %>% 
  group_by(tissue, time, sample, clonotype_id_tcr) %>%   
  tally() %>% 
  rename(freq = "n")

# 1) calculate the Gini coefficient per cell_type & per sample
gini_res <- clons_for_gini %>% 
  group_by(tissue, time, sample) %>%        ## tissue, time,
  summarise(across(-clonotype_id_tcr, Gini, na.rm = TRUE)) %>% 
  filter(!freq %in% NaN)

gini_res <- gini_res %>% 
  filter(!(tissue %in% 'CSF' & time %in% 'YR2'))

# 1) order the labels
# 2) specify the consistent colors
gini_res$time = factor(gini_res$time, levels=c('BSE', 'W5_10', 'YR1', 'YR2'))    
gini_res$tissue = factor(gini_res$tissue, levels=c('PBMC', 'CSF'))

# 1) create a column "subject"
gini_res$sample <- as.character(gini_res$sample)
gini_res['subject'] <- sapply(gini_res$sample, function(x){
  strsplit(x, split = "_", fixed = T)[[1]][1]})

# 1) run tukey test (NORMAL)
tukey_stat <- gini_res %>%
  group_by(tissue) %>%
  tukey_hsd(freq ~ time, p.adjust.method = 'fdr') %>%     # tissue; time
  filter(group1 == 'BSE') %>%          # PBMC; BSE
  add_xy_position(x = "time") 

# 1) specify your pallete
# BSE: #ea698b, W5_10: #dec9e9, YR1: #a06cd5, YR2: #47126b
my_palette <- c('#ea698b', '#dec9e9', '#a06cd5', '#47126b')    

# 1) BOXPLOTS WITH STATS: plot the values of Gini coefficient per cell type
fig <- ggplot(gini_res, aes(x = time, y = freq)) +    #fig <- 
  geom_boxplot(outlier.shape = NA, aes(fill = time), alpha = 0.7) +    # tissue; alpha = 0.5
  geom_point(position=position_dodge(width=0.75), alpha = 0.7,        # alpha = 0.5
             aes(fill = time), colour="black", pch=21) +             # tissue
  geom_line(aes(x = as.numeric(factor(time)) + .75 / 4 * ifelse(time == "BSE", -1, 1),   # tissue == "PBMC"; time == "BSE"
                group = interaction(subject, time)), color = "black", linewidth = 0.5, alpha = 0.5) +
  ylab("Gini coefficient") +
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
  # scale_fill_discrete(name = "Tissue") +
  scale_fill_manual(name = "Time", values = my_palette) +
  ggtitle(paste0("TCR repertoire diversity of CD8 T cells across timepoints")) +
  facet_grid(~tissue)

# 1) save in this to keep the good asterisk shape
ggsave(filename=paste0('output/figures/TCR/Gini_coef/CD4_CD8/CD8_PBMC_CSF.pdf'), 
       width = 6, height = 3, dpi = 300)
png(filename=paste0('output/figures/TCR/Gini_coef/CD4_CD8/CD8_PBMC_CSF.png'), 
    width = 6, height = 3, units = "in", res = 300)   # , org, '_BSE_', time_filt, '.png'
plot(fig)
dev.off()

#-------------------------------Antigen prediction based on CDR3_aa--------------------

# 1) load the antigen database
vdj_db <- read.csv(file = 'data/TCR/VDJ_db.tsv', sep = '\t')

# 1) remove useless columns
vdj_db_s <- vdj_db %>% 
  select(Gene, CDR3, Species, MHC.A, MHC.B, Epitope, Epitope.gene, Epitope.species) %>% 
  rename(Chain = Gene, CDR3_aa = CDR3) %>% 
  mutate(chain.cdr3_aa = paste0(Chain, ':', CDR3_aa)) %>% 
  filter(Species %in% 'HomoSapiens')

# 1) split vdj_db_s into TRA & TRB based seqs
vdj_alpha <- vdj_db_s %>% filter(Chain %in% "TRA") %>% rename(TRA_chain = chain.cdr3_aa)
vdj_betta <- vdj_db_s %>% filter(Chain %in% "TRB") %>% rename(TRB_chain = chain.cdr3_aa)

# 1) split "chain.cdr3_aa" into TRA & TRB chain
# 2) swap TRA sequnces in TRB_chain column and visa versa
trc_df <- expData@meta.data %>% 
  separate(cdr3s_aa, c('TRA_chain', "TRB_chain"), sep=";", remove = F)
trc_df[!(substr(trc_df$TRA_chain, 1, 3) %in% "TRA"), c("TRB_chain", "TRA_chain") ] <- trc_df[!(substr(trc_df$TRA_chain, 1, 3) %in% "TRA"), c("TRA_chain", "TRB_chain") ]

# 1) map vdj_betta database to GEX-TCR meta
mapped_b <- merge(trc_df, vdj_betta, by = "TRB_chain", all.x=T)
mapped_b <- mapped_b[!duplicated(mapped_b$barcode), ]

# 1) Leave cells with TCR=True & Chain =NA for vdj_alpha mapping
unmapped_a <- mapped_b %>% 
  filter(contains_tcr %in% TRUE & is.na(Chain)) %>% 
  select(-c('Chain', 'CDR3_aa', 'Species', 'MHC.A',
            'MHC.B', 'Epitope', 'Epitope.gene', 'Epitope.species'))

# 1) map vdj_alpha database to GEX-TCR meta
mapped_a <- merge(unmapped_a, vdj_alpha, by = "TRA_chain", all.x=T)
mapped_a <- mapped_a[!duplicated(mapped_a$barcode), ]

# 1) Remove cells with TCR=True & Chain =NA from vdj_betta df
mapped_b <- mapped_b %>% 
  filter(!(contains_tcr %in% TRUE & is.na(Chain)))

# 1) merge two DFs together by rows
# 2) restore the order of cells
mapped_vdj_db <- rbind(mapped_b, mapped_a)
rownames(mapped_vdj_db) <- mapped_vdj_db$barcode
mapped_vdj_db <- mapped_vdj_db[order(mapped_vdj_db$ordered_id), ]

# 1) save and reload
path_3 <- '/storage1/fs1/martyomov/Active/collaborations/romans/Greg_wu/CLOCK/output/output/TCR/'
write.csv(mapped_vdj_db, paste0(path_3, 'mapped_vdj_db.csv'), row.names=T)
mapped_vdj_db <- read.csv(paste0(path_3, 'mapped_vdj_db.csv'), row.names = 1)

# 1) count the #cells with TCRs grouped by epitope species
ag_freq <- mapped_vdj_db %>% 
  filter(contains_tcr %in% TRUE) %>% 
  group_by(Epitope.species) %>% 
  tally() %>% 
  rename(freq = n) %>% 
  mutate(Epitope.species = replace_na(Epitope.species, "Unknown")) %>% 
  filter(!Epitope.species %in% c("Homo sapiens", "HomoSapiens", "M.tuberculosis", "PseudomonasFluorescens", "Wheat"))

# 1) plot 
ggplot(ag_freq, aes(x = reorder(Epitope.species, -freq), y = freq)) +
  geom_bar(stat="identity") +
  theme_light() +
  xlab("Epitope species") +
  ylab("# clones") +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  scale_y_break(c(3222, 31800))     # can add 2nd break: scale_y_break(c(763, 3000))
path <- '/storage1/fs1/martyomov/Active/collaborations/romans/Greg_wu/CLOCK/output/'
ggsave(paste0(path, 'figures/TCR/Freq_barplot/antigen_pred/all_epit_sp.png'), width = 10, height = 7)
ggsave(paste0(path, 'figures/TCR/Freq_barplot/antigen_pred/all_epit_sp.pdf'), width = 10, height = 7)

# 1) count the #cells with TCRs grouped by epitope species & tissue
ag_freq <- mapped_vdj_db %>% 
  filter(contains_tcr %in% TRUE) %>% 
  group_by(Epitope.species, Chain) %>%      #Chain/tissue
  tally() %>% 
  rename(freq = n) %>% 
  mutate(Epitope.species = replace_na(Epitope.species, "Unknown")) %>% 
  filter(!Epitope.species %in% c("Homo sapiens", "HomoSapiens", "M.tuberculosis", "PseudomonasFluorescens", "Wheat"))

# 1) order the labels
# ag_freq$tissue = factor(ag_freq$tissue, levels=c('PBMC', 'CSF'))
ag_freq$Chain = factor(ag_freq$Chain, levels=c('TRA', 'TRB'))

# 1) plot 
axis_text_size = 16
title_text_size = 18
ggplot(ag_freq, aes(x = reorder(Epitope.species, -freq), y = freq)) +
  geom_bar(stat="identity", aes(fill = Chain), alpha = 0.8) +        #Chain/tissue
  theme_light() +
  theme(axis.text.x = element_text(size = axis_text_size),
        axis.text.y = element_text(size = axis_text_size),
        axis.title.x = element_text(size = title_text_size),
        axis.title.y = element_text(size = title_text_size),
        legend.text = element_text(size = axis_text_size),
        legend.title = element_text(size = title_text_size)) +
  xlab("Viral epitops") +
  ylab("# clones") +
  scale_fill_discrete(name = "Chain") +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  scale_y_break(c(3222, 31800))     # can add 2nd break: scale_y_break(c(763, 3000))
ggsave(paste0(path, 'figures/TCR/Freq_barplot/antigen_pred/all_epit_sp_chain.png'), 
       width = 10, height = 7)
ggsave(paste0(path, 'figures/TCR/Freq_barplot/antigen_pred/all_epit_sp_chain.pdf'), 
       width = 10, height = 7)

#-------------------------------------Analysis of viral-specific subsets-----------------

# 1) subset only virus-specific clones
viruses <- mapped_vdj_db %>% 
  filter(contains_tcr %in% TRUE & Epitope.species %in% c("SARS-CoV-2", "CMV", "EBV", "HCV", "InfluenzaA", "HIV-1"                 
                                                         ,"YFV" ,"HPV", "HCoV-HKU1", "HTLV-1", "DENV3/4", "HSV-2"))

### 1) remove useless columns FOR SIMPLE BOXPLOT
# 2) and count clone sizes
org <- "PBMC"
clons_for_gini <- viruses %>% 
  filter(tissue %in% org &
           time %in% "BSE" &
           cell_type %in% c('CD4 Naive', 'CD4 Treg naive', 'CD4 Treg memory', 'CD4 Th1', 'CD4 Th17', 'CD4 Th22', 'CD4 Temra', 'CD4 exhausted',
                            'CD8 Naive', 'CD8 Tcm CCR4-', 'CD8 Tcm CCR4+', 'CD8 Trm', 'CD8 Tem GZMK+', 'CD8 Tem GZMB+', 'CD8 Temra', 'CD8 NKT-like') ) %>%   
  filter(contains_tcr == TRUE) %>% 
  group_by(cell_type, sample, clonotype_id) %>%    
  tally() %>% 
  rename(freq = "n")

# 1) calculate the Gini coefficient per cell_type & per sample
gini_res <- clons_for_gini %>% 
  group_by(cell_type, sample) %>%       
  summarise(across(-clonotype_id, Gini, na.rm = TRUE)) %>% 
  filter(!freq %in% NaN)

# 1) order the labels
# 2) specify the consistent colors
gini_res$cell_type = factor(gini_res$cell_type, levels=c('CD4 Naive', 'CD4 Treg naive', 'CD4 Treg memory', 'CD4 Th1', 'CD4 Th17', 'CD4 Th22', 'CD4 Temra', 'CD4 exhausted',
                                                         'CD8 Naive', 'CD8 Tcm CCR4-', 'CD8 Tcm CCR4+', 'CD8 Trm', 'CD8 Tem GZMK+', 'CD8 Tem GZMB+', 'CD8 Temra', 'CD8 NKT-like'))
my_palette <- c('#d60000',  '#ff7266', '#edb8b8', '#e36414', '#ff9f1c', '#ffbf69', '#f9cb9cff', '#fdf490', 
                '#a3b18a', '#588157', '#3a5a40', '#bce784', '#55a630', '#007f5f', '#a5be00', '#74d3ae')

## 1) BOXPLOTS WITHOUT STATS: plot the values of Gini coefficient per cell type
ggboxplot(gini_res, x = 'cell_type', y = 'freq', fill = 'cell_type') +
  geom_point(position=position_dodge(width=0.8),
             aes(fill = cell_type), colour="black", pch=21) +
  ylab(paste0("Gini coefficient")) +
  theme_light() +
  theme(axis.title.x = element_blank(),
        axis.title.y=element_text(size=10),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        legend.position="none",
        plot.title = element_text(size = 10, hjust = 0.5))   +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  scale_fill_manual(values = my_palette) +
  ggtitle(paste0("TCR repertoire diversity in ", org))

ggsave(paste0('output/figures/TCR/Gini_coef/virus_specific/CD4_CD8_', org, '_all.png'),
       width = 7, height = 4,  dpi = 300)
ggsave(paste0('output/figures/TCR/Gini_coef/virus_specific/CD4_CD8_', org, '_all.pdf'),
       width = 7, height = 4,  dpi = 300)

### 1) remove useless columns FOR COMPLEX/ADVANCED BOXPLOT
# 2) and count clone sizes
org <- c('PBMC', 'CSF')
# time_filt <- 'YR2'
clons_for_gini <- viruses %>% 
  mutate(time = recode(time, 'W05' = 'W5_10', 'W10' = 'W5_10')) %>%   #'YR1' = 'YR1_2', 'YR2' = 'YR1_2',
  filter(tissue %in% org &
           cell_type %in% c('CD4 Naive', 'CD4 Treg naive', 'CD4 Treg memory', 'CD4 Th1', 'CD4 Th17', 'CD4 Th22', 'CD4 Temra', 'CD4 exhausted',
                            'CD8 Naive', 'CD8 Tcm CCR4-', 'CD8 Tcm CCR4+', 'CD8 Trm', 'CD8 Tem GZMK+', 'CD8 Tem GZMB+', 'CD8 Temra', 'CD8 NKT-like') &
           time %in% c("BSE")                                                ##
  ) %>%   
  group_by(cell_type, tissue, sample, clonotype_id) %>%    # tissue, time
  tally() %>% 
  rename(freq = "n")

# 1) calculate the Gini coefficient per cell_type & per sample
gini_res <- clons_for_gini %>% 
  group_by(cell_type, tissue, sample) %>%        ## tissue, time
  summarise(across(-clonotype_id, Gini, na.rm = TRUE)) %>% 
  filter(!freq %in% NaN)

# 1) order the labels
# 2) specify the consistent colors
gini_res$cell_type = factor(gini_res$cell_type, levels=c('CD4 Naive', 'CD4 Treg naive', 'CD4 Treg memory', 'CD4 Th1', 'CD4 Th17', 'CD4 Th22', 'CD4 Temra', 'CD4 exhausted',
                                                         'CD8 Naive', 'CD8 Tcm CCR4-', 'CD8 Tcm CCR4+', 'CD8 Trm', 'CD8 Tem GZMK+', 'CD8 Tem GZMB+', 'CD8 Temra', 'CD8 NKT-like'))
# gini_res$time = factor(gini_res$time, levels=c('BSE', time_filt))
gini_res$tissue = factor(gini_res$tissue, levels=c('PBMC', 'CSF'))

# 1) create a column "subject"
gini_res$sample <- as.character(gini_res$sample)
gini_res['subject'] <- sapply(gini_res$sample, function(x){
  strsplit(x, split = "_", fixed = T)[[1]][1]})

# 1) run tukey test (NORMAL)
tukey_stat <- gini_res %>%
  filter(!cell_type %in% "CD8 Trm") %>% 
  group_by(cell_type) %>%
  tukey_hsd(freq ~ tissue, p.adjust.method = 'fdr') %>%     # tissue, time
  filter(group1 == 'PBMC') %>%          # PBMC; BSE
  add_xy_position(x = "cell_type") 

tukey_stat <- tukey_stat %>% 
  filter(!cell_type %in% c("CD4 exhausted", 'CD8 NKT-like'))
gini_res <- gini_res %>% 
  filter(!cell_type %in% c("CD4 exhausted", 'CD8 NKT-like'))

# 1) specify your pallete
# BSE: #ea698b, W5_10: #dec9e9, YR1: #a06cd5, YR2: #47126b
# my_palette <- c('#ea698b', '#47126b')    

# 1) BOXPLOTS WITH STATS: plot the values of Gini coefficient per cell type
fig <- ggplot(gini_res, aes(x = cell_type, y = freq)) +
  geom_boxplot(outlier.shape = NA, aes(fill = tissue), alpha = 0.7) +    # tissue; alpha = 0.5
  geom_point(position=position_dodge(width=0.75), alpha = 0.7,        # alpha = 0.5
             aes(fill = tissue), colour="black", pch=21) +             # tissue
  geom_line(aes(x = as.numeric(factor(cell_type)) + .75 / 4 * ifelse(tissue == "PBMC", -1, 1),   # tissue == "PBMC"
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
  scale_fill_discrete(name = "Tissue") +     #, values = my_palette
  ggtitle(paste0("Virus-specific TCR repertoire diversity in BSE: PBMC vs CSF"))

# 1) save in this to keep the good asterisk shape
ggsave('output/figures/TCR/Gini_coef/virus_specific/PBMC_CSF_BSE.pdf', width = 10, height = 4, dpi = 300)
png(filename=paste0('output/figures/TCR/Gini_coef/virus_specific/PBMC_CSF_BSE.png'), 
    width = 10, height = 4, units = "in", res = 300)   # open the device
plot(fig)
dev.off()

#-------------------------------------DONUT CHART: CLADT & VIRAL CLONALITY-----------------

# 1) subset only virus-specific clones
viruses <- mapped_vdj_db %>% 
  filter(contains_tcr %in% TRUE & Epitope.species %in% c("SARS-CoV-2", "CMV", "EBV", "HCV", "InfluenzaA", "HIV-1"                 
                                                         ,"YFV" ,"HPV", "HCoV-HKU1", "HTLV-1", "DENV3/4", "HSV-2"))

# 1) rename treatment groups into: Before CladT vs After CladT
viruses <- viruses %>% 
  mutate(time = recode(time, 'BSE'='Before CladT',
                       'W05'= 'After CladT', 'W10'='After CladT', 
                       'YR1'='After CladT', 'YR2'='After CladT'))

# 1) DONUT CHART: calculate frequncies of clonotypes grouped by tissue & cell_type_2 & pick the biggest
freq_tcr <- viruses %>%
  mutate(time_comb = time) %>% 
  group_by(cell_type_2, tissue, time_comb, clonotype_id) %>%    #cell_type
  tally()  

# 1) keep only those groups in which you are interested in
# organ <- 'PBMC'
# time <- 'YR2'
freq_tcr_filt <- freq_tcr %>% 
  filter(tissue %in% 'CSF' & 
         cell_type_2 %in% 'CD4 T cell') %>%   #cell_type_2 %in% 'CD8 T cell'
  rename(clon_size = "n")
freq_tcr_filt$clon_size %>% table()

# 1) count the total #clones per time group
# 2) calculate in % the contribution of clones in each group
num_tot_cells <- freq_tcr_filt %>% 
  group_by(time_comb) %>%   #tissue, time_comb
  summarise(total_clon_size = sum(clon_size)) 
freq_tcr_filt$total_clon_size <- with(freq_tcr_filt, 
                                      num_tot_cells$total_clon_size[match(time_comb, num_tot_cells$time_comb)])   #time_comb, tissue
freq_tcr_filt <- freq_tcr_filt %>% 
  mutate(percent = clon_size/total_clon_size * 100)

# 1) create the colomn "clone_group" using when_case & split into 4 groups (unique, 2-4, 5-9, >=10)
freq_tcr_filt <- freq_tcr_filt %>% 
  mutate(clone_group = case_when(clon_size == 1 ~ 'Unique',        ## for Greg: clonotype_id == 'clonotype_24487' ~ 'clonotype_24487',
                                 clon_size >=2 & clon_size <= 4 ~ '2-4',
                                 clon_size >=5 & clon_size <= 9 ~ '5-9',
                                 clon_size >= 10 ~ '>=10'))

# 1) collapse based on clone_group/Unique group
freq_tcr_filt <- freq_tcr_filt %>% 
  group_by(time_comb,    #cell_type_2, tissue,
           grp = case_when(clone_group %in% "Unique" ~ 0, TRUE ~ row_number())) %>% 
  summarise(across(starts_with(c("clon_size", "percent")), sum),
            clone_group = str_c(clone_group, collapse=",")) %>% 
  select(-grp) %>% 
  mutate(clone_group = gsub("^Uni.*", "Unique", clone_group))

# 1) collapse based on "2-4" group 
freq_tcr_filt <- freq_tcr_filt %>% 
  group_by(time_comb,    #cell_type_2, tissue,
           grp = case_when(clone_group %in% "2-4" ~ 0, TRUE ~ row_number())) %>% 
  summarise(across(starts_with(c("clon_size", "percent")), sum),
            clone_group = str_c(clone_group, collapse=",")) %>% 
  select(-grp) %>% 
  mutate(clone_group = gsub("^2.*", "2-4", clone_group))

# 1) order the labels
freq_tcr_filt$clone_group %>% unique()
freq_tcr_filt$clone_group = factor(freq_tcr_filt$clone_group,
                                   levels=c("Unique", "2-4", "5-9", ">=10"))   #, "clonotype_24487"
freq_tcr_filt$time_comb = factor(freq_tcr_filt$time_comb, levels=c("Before CladT", "After CladT"))
# freq_tcr_filt$tissue = factor(freq_tcr_filt$tissue, levels=c("PBMC", "CSF"))

# 1) plot the Donut chart
colors <- c("#adb5bd", "#ff8fa3", "#ff4d6d", "#a4133c")   #, "#e09f3e"
ggplot(data = freq_tcr_filt, aes(x=2, y = percent, fill = clone_group))+    #, fill = clone_group
  # geom_bar(stat="identity", color='black')+
  # geom_bar(stat="identity", aes(fill = clone_group))+
  geom_col(color='black')+
  coord_polar("y", start = 0) + 
  facet_wrap(~ time_comb) +    #tissue, time_comb
  theme(panel.background = element_blank(),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(), 
        strip.text = element_text(size = 15)) +
  scale_fill_manual(values = colors, name = 'No. of clones') +
  xlim(0.5, 2.5)  
ggsave('output/figures/TCR/Donut_chart/across_all/virus/CD4_CSF_cladt.png', width = 10, height = 7)
ggsave('output/figures/TCR/Donut_chart/across_all/virus/CD4_CSF_cladt.pdf', width = 10, height = 7)


#----------------Dynamics of T cell migration (aka Brian's experiment): NOT NORMILIZED----------------------

# 1) set env & load the data
path <- "/storage1/fs1/martyomov/Active/collaborations/romans/Greg_wu/CLOCK/"
setwd(path)
expData <- readRDS('output/output/scRNA/combined_annot.rds')

# 1) subset the df & select only useful comns
tcr_db <- expData@meta.data %>% 
  filter(contains_tcr %in% TRUE) %>% 
  mutate(time = recode(time, 'W05' = 'W5/10', 'W10' = 'W5/10')) %>% 
  select(subject, time, sample, tissue, clonotype_id_tcr) %>%    #barcode,
  filter(sample %in% c('005_BSE_CSF', '005_BSE_PBMC', '006_BSE_CSF', '006_BSE_PBMC', '009_BSE_CSF', '009_BSE_PBMC'
                       , '013_BSE_CSF', '013_BSE_PBMC'                 # '011_BSE_CSF', '011_BSE_PBMC',
                       , '0015_BSE_CSF', '0015_BSE_PBMC', '011_W10_CSF', '011_W10_PBMC'
                       , '013_W05_CSF', '013_W05_PBMC', '005_YR1_CSF', '005_YR1_PBMC'
                       , '006_YR1_CSF', '006_YR1_PBMC'))   #, '001_YR2_CSF', '001_YR2_PBMC'

# 1) remove singletons (clone size =1)
tcr_nsingl <- subset(tcr_db, duplicated(clonotype_id_tcr) | duplicated(clonotype_id_tcr,
                                                                       fromLast=TRUE))

# 1) group clonotypes by subj, clonotype_id, time and tissue & count freq
tcr_grp <- tcr_nsingl %>% 
  group_by(subject, clonotype_id_tcr, time, tissue) %>% 
  tally() # %>% filter(n > 1)

# 1) count their clone_sizes within a subj & time only
clons_size <- tcr_nsingl %>% 
  group_by(subject, time, clonotype_id_tcr) %>% 
  tally() %>% 
  rename(clon_size = n)

# 1) count total # unique clonotypes per subj, per time 
tcr_nsc_tot <- tcr_grp %>% 
  group_by(subject, time) %>% 
  distinct(clonotype_id_tcr) 

# 1) define their size by merging with "clons_size" DF
# 2) filter clones that are expanded in a subject within a timepoint
tcr_nsc_tot_lst <- merge(tcr_nsc_tot, clons_size, 
                         by=c("subject", "time", "clonotype_id_tcr"), all.x=TRUE)
tcr_nsc_tot_lst <- tcr_nsc_tot_lst %>% 
  filter(clon_size > 1)

# 1) sum clone size values
tcr_nsc_tot_s <- tcr_nsc_tot_lst %>% 
  group_by(subject, time) %>% 
  tally(clon_size) %>% 
  rename(NSC_tot = n) 

# 1) leave only unique clonotypes per subj, per time that are in PBMC & CSF
tcr_nsc_p.c_lst <- tcr_grp %>% 
  group_by(subject, time, clonotype_id_tcr) %>% 
  tally() %>% 
  filter(n == 2) %>% 
  group_by(subject, time) %>% 
  distinct(clonotype_id_tcr) 

# 1) define their size by merging with "clons_size" DF
# 2) sum clone size values
tcr_nsc_p.c_lst <- merge(tcr_nsc_p.c_lst, clons_size, by=c("subject", "time", "clonotype_id_tcr"), all.x=TRUE)
tcr_nsc_p.c_s <- tcr_nsc_p.c_lst %>% 
  group_by(subject, time) %>% 
  tally(clon_size) %>%        # similar to "sum(clon_size)"
  rename(NSC_p.c = n)

# 1) merge two DFs
# 2) count the precent per subject & time
tcr_df_m <- merge(tcr_nsc_p.c_s, tcr_nsc_tot_s, by=c("subject", "time"), all.x=TRUE)
tcr_df_m <- tcr_df_m %>% mutate(percent = NSC_p.c/NSC_tot * 100)

# 1) run tukey test (NORMAL)
tukey_stat <- tcr_df_m %>%
  tukey_hsd(percent ~ time, p.adjust.method = 'fdr') %>%     # tissue; time
  filter(group1 == 'BSE') %>%         
  add_xy_position(x = "time") 

# 1) specify your palette
time_palette <- c('#ea698b', '#dec9e9', '#a06cd5')      #BSE, W5_10, YR1

# 1) calculate mean and sd
tcr_df_grp <- tcr_df_m %>% 
  group_by(time) %>% 
  summarise(mean=mean(percent), sd=sd(percent))

# 1) plot line graph
ggplot(tcr_df_grp, aes(x = time, y = mean)) +
  geom_line(aes(group=1, color=time)) +
  geom_point(aes(color=time)) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd, color = time), width=.2,
                position=position_dodge(0.05)) +
  expand_limits(y = 0) +
  ylab("% non-singleton clones in both PBMC & CSF") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y=element_text(size=10),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        legend.position="right",
        legend.title=element_text(size=9),
        legend.text=element_text(size=8),
        plot.title = element_text(size = 10, hjust = 0.5))+
  stat_pvalue_manual(tukey_stat, hide.ns = T) +
  scale_color_manual(name = "Time", values = time_palette) +
  ggtitle(paste0("Dynamics of T cell migration after cladribine treatment"))

# save
ggsave(filename=paste0('output/figures/TCR/Line_plot/migration_norm_strict.pdf'), 
       width = 5, height = 3.5, dpi = 300)
ggsave(filename=paste0('output/figures/TCR/Line_plot/migration_norm_strict.png'), 
       width = 6, height = 4, dpi = 300)

# 1) BOXPLOTS WITH STATS: 
fig <- ggplot(tcr_df_m, aes(x = time, y = percent)) +  #fig <- 
  geom_boxplot(outlier.shape = NA, aes(fill = time), alpha = 0.7) +    # tissue; alpha = 0.5
  geom_point(position=position_dodge(width=0.75), alpha = 0.7,        # alpha = 0.5
             aes(fill = time), colour="black", pch=21) +          
  ylab("% non-singleton clones in both PBMC & CSF") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y=element_text(size=10),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        legend.position="right",
        legend.title=element_text(size=9),
        legend.text=element_text(size=8),
        plot.title = element_text(size = 10, hjust = 0.5))+
  stat_pvalue_manual(tukey_stat, hide.ns = T) +
  scale_fill_manual(name = "Time", values = time_palette) +
  ggtitle(paste0("Dynamics of T cell migration after cladribine treatment"))   #"TCR repertoire diversity in BSE: PBMC vs CSF"

# 1) save in this to keep the good asterisk shape
ggsave(filename=paste0('output/figures/TCR/Freq_barplot/migration/unnorm_no_YR2_no_Singl.pdf'), 
       width = 5, height = 3.5, dpi = 300)
png(filename=paste0('output/figures/TCR/Freq_barplot/migration/unnorm_no_YR2_no_Singl.png'), 
    width = 5, height = 3.5, units = "in", res = 300)  
plot(fig)
dev.off()

#----------------Dynamics of T cell migration: NORMILIZED BY CLONE SIZE----------------------

# 1) set env & load the data
path <- "/storage1/fs1/martyomov/Active/collaborations/romans/Greg_wu/CLOCK/"
setwd(path)
expData <- readRDS('output/output/scRNA/combined_annot.rds')

# 1) subset the df & select only useful comns
tcr_db <- expData@meta.data %>% 
  filter(contains_tcr %in% TRUE) %>% 
  mutate(time = recode(time, 'W05' = 'W5/10', 'W10' = 'W5/10')) %>% 
  select(subject, time, sample, tissue, clonotype_id_tcr) %>%    #barcode,
  filter(sample %in% c('005_BSE_CSF', '005_BSE_PBMC', '006_BSE_CSF', '006_BSE_PBMC', '009_BSE_CSF', '009_BSE_PBMC'
                       , '013_BSE_CSF', '013_BSE_PBMC'                 # '011_BSE_CSF', '011_BSE_PBMC',
                       , '0015_BSE_CSF', '0015_BSE_PBMC', '011_W10_CSF', '011_W10_PBMC'
                       , '013_W05_CSF', '013_W05_PBMC', '005_YR1_CSF', '005_YR1_PBMC'
                       , '006_YR1_CSF', '006_YR1_PBMC'))   #, '001_YR2_CSF', '001_YR2_PBMC'

# 1) remove singletons (clone size =1)
tcr_nsingl <- subset(tcr_db, duplicated(clonotype_id_tcr) | duplicated(clonotype_id_tcr,
                                                                       fromLast=TRUE))

# 1) group clonotypes by subj, clonotype_id, time and tissue & count freq
tcr_grp <- tcr_nsingl %>% 
  group_by(subject, clonotype_id_tcr, time, tissue) %>% 
  tally() # %>% filter(n > 1)

# 1) count their clone_sizes within a subj & time only
clons_size <- tcr_nsingl %>% 
  group_by(subject, time, clonotype_id_tcr) %>% 
  tally() %>% 
  rename(clon_size = n)

# 1) count total # unique clonotypes per subj, per time 
tcr_nsc_tot <- tcr_grp %>% 
  group_by(subject, time) %>% 
  distinct(clonotype_id_tcr) 

# 1) define their size by merging with "clons_size" DF
# 2) sum clone size values
tcr_nsc_tot_lst <- merge(tcr_nsc_tot, clons_size, by=c("subject", "time", "clonotype_id_tcr"), all.x=TRUE)
tcr_nsc_tot_s <- tcr_nsc_tot_lst %>% 
  group_by(subject, time) %>% 
  tally(clon_size) %>% 
  rename(NSC_tot = n)

# 1) leave only unique clonotypes per subj, per time that are in PBMC & CSF
tcr_nsc_p.c_lst <- tcr_grp %>% 
  group_by(subject, time, clonotype_id_tcr) %>% 
  tally() %>% 
  filter(n == 2) %>% 
  group_by(subject, time) %>% 
  distinct(clonotype_id_tcr) 

# 1) define their size by merging with "clons_size" DF
# 2) sum clone size values
tcr_nsc_p.c_lst <- merge(tcr_nsc_p.c_lst, clons_size, by=c("subject", "time", "clonotype_id_tcr"), all.x=TRUE)
tcr_nsc_p.c_s <- tcr_nsc_p.c_lst %>% 
  group_by(subject, time) %>% 
  tally(clon_size) %>%        # similar to "sum(clon_size)"
  rename(NSC_p.c = n)
  
# 1) merge two DFs
# 2) count the precent per subject & time
tcr_df_m <- merge(tcr_nsc_p.c_s, tcr_nsc_tot_s, by=c("subject", "time"), all.x=TRUE)
tcr_df_m <- tcr_df_m %>% mutate(percent = NSC_p.c/NSC_tot * 100)

# 1) run tukey test (NORMAL)
tukey_stat <- tcr_df_m %>%
  tukey_hsd(percent ~ time, p.adjust.method = 'fdr') %>%     # tissue; time
  filter(group1 == 'BSE') %>%         
  add_xy_position(x = "time") 

# 1) run dunn's test (NOT NORMAL)
dunn_stat <- tcr_df_m %>%
  dunn_test(percent ~ time, p.adjust.method = 'fdr') %>%
  filter(group1 == 'BSE') %>%
  add_xy_position(x = "time")

# 1) specify your palette
time_palette <- c('#ea698b', '#dec9e9', '#a06cd5')      #BSE, W5_10, YR1

# 1) calculate mean and sd
tcr_df_grp <- tcr_df_m %>% 
  group_by(time) %>% 
  summarise(mean=mean(percent), sd=sd(percent))

# 1) plot line graph
ggplot(tcr_df_grp, aes(x = time, y = mean)) +
  geom_line(aes(group=1, color=time)) +
  geom_point(aes(color=time)) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd, color = time), width=.2,
                position=position_dodge(0.05)) +
  expand_limits(y = 0) +
  ylab("% non-singleton clones in both PBMC & CSF") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y=element_text(size=10),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        legend.position="right",
        legend.title=element_text(size=9),
        legend.text=element_text(size=8),
        plot.title = element_text(size = 10, hjust = 0.5))+
  stat_pvalue_manual(tukey_stat, hide.ns = T) +
  scale_color_manual(name = "Time", values = time_palette) +
  ggtitle(paste0("Dynamics of T cell migration after cladribine treatment"))

# save
ggsave(filename=paste0('output/figures/TCR/Line_plot/migration_norm.pdf'), 
       width = 5, height = 3.5, dpi = 300)
ggsave(filename=paste0('output/figures/TCR/Line_plot/migration_norm.png'), 
       width = 6, height = 4, dpi = 300)

# 1) BOXPLOTS WITH STATS: 
ggplot(tcr_df_m, aes(x = time, y = percent)) +  #fig <- 
  geom_boxplot(outlier.shape = NA, aes(fill = time), alpha = 0.7) +    # tissue; alpha = 0.5
  geom_point(position=position_dodge(width=0.75), alpha = 0.7,        # alpha = 0.5
             aes(fill = time), colour="black", pch=21) +          
  ylab("% non-singleton clones in both PBMC & CSF") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y=element_text(size=10),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        legend.position="right",
        legend.title=element_text(size=9),
        legend.text=element_text(size=8),
        plot.title = element_text(size = 10, hjust = 0.5))+
  stat_pvalue_manual(tukey_stat, hide.ns = T) +
  scale_fill_manual(name = "Time", values = time_palette) +
  ggtitle(paste0("Dynamics of T cell migration after cladribine treatment"))   #"TCR repertoire diversity in BSE: PBMC vs CSF"

# 1) save in this to keep the good asterisk shape
ggsave(filename=paste0('output/figures/TCR/Freq_barplot/migration/norm_no_YR2_no_Singl.pdf'), 
       width = 5, height = 3.5, dpi = 300)
png(filename=paste0('output/figures/TCR/Freq_barplot/migration/norm_no_YR2_no_Singl.png'), 
    width = 5, height = 3.5, units = "in", res = 300)  
plot(fig)
dev.off()

#-----------------------------Follow-ups on Brian's comments------------------

# 1) set env & load the data
path <- "/storage1/fs1/martyomov/Active/collaborations/romans/Greg_wu/CLOCK/"
setwd(path)
expData <- readRDS('output/output/scRNA/combined_annot.rds')
barcodes_f <- read.csv("output/output/TCR/filtered_TCR.csv")
contigs_unf <- read.csv("output/output/TCR/unfiltered_TCR.csv")

# 1) standardize barcode names
contigs_unf <- contigs_unf %>% 
  mutate(barcode = gsub("-1", "", unique_barcode))

# 1) leave only cells associated with "clonotype_7449"
barcodes_f_s <- barcodes_f %>% 
  filter(clonotype_id %in% 'clonotype_7449')

# 1) subset by barcodes (should be 362 rows)
clon_7449_tcr <- contigs_unf %>% 
  filter((barcode %in% barcodes_f_s$barcode) & !(productive %in% 'False')) %>% 
  select(barcode, contig_id, chain, v_gene, j_gene, cdr3, cdr3_nt)

# 1) compress it by barcode & add chain, V and J gene colmns (NO CDR3_NT)
clon_7449_tcr <- clon_7449_tcr %>% 
  mutate(chain.cdr3_nt = paste0(chain, ':', cdr3_nt)) %>% 
  mutate(chain_v_gene = paste0(chain, ':', v_gene)) %>% 
  mutate(chain_j_gene = paste0(chain, ':', j_gene)) %>% 
  select(barcode, chain.cdr3_nt, chain_v_gene, chain_j_gene) %>%  
  group_by(barcode) %>% 
  mutate(cdr3s_nt = paste0(chain.cdr3_nt, collapse = '_')) %>% 
  mutate(v_genes = paste0(chain_v_gene, collapse = '_')) %>% 
  mutate(j_genes = paste0(chain_j_gene, collapse = '_')) %>% 
  mutate(cdr3nt_v_j = paste0(cdr3s_nt, v_genes, j_genes, collapse = '_')) %>% 
  select(-c(chain.cdr3_nt, chain_v_gene, chain_j_gene)) %>% 
  unique() 

# 1) redifine clonotypes based on cdr3aa, v_gene & j_gene
clon_7449_tcr <- clon_7449_tcr %>% 
  group_by(cdr3nt_v_j) %>% 
  mutate(clonotype_id = cur_group_id()) %>% 
  mutate(clonotype_id = paste0('clonotype_', clonotype_id))

# 1) leave only cells associated with "clonotype_7449"
# 2) remove useless colmns
clon_7449_gex <- expData@meta.data %>% 
  filter(clonotype_id_tcr %in% 'clonotype_7449') %>% 
  select(barcode, cell_type)

# 1) map "cell_type" colmn to clon_7449_tcr df
clon_7449_cell_type <- merge(clon_7449_tcr, clon_7449_gex, by=c("barcode"), all.x=TRUE)
write.csv(clon_7449_cell_type, file='output/output/TCR/cd4_cd8_clon/cdr3nt_v_j.csv')

# 1) count #clones per tissue, time & cell_type
clon_7449_freq <- clon_7449_cell_type %>% 
  group_by(clonotype_id, cell_type) %>% 
  tally() %>% 
  rename(clone_size = n)

# 1) fix the order of labels
clon_7449_freq$cell_type <- factor(clon_7449_freq$cell_type,
                                      levels=c('CD4 Naive', 'CD4 Treg naive', 'CD4 Treg memory', 'CD4 Th1',
                                               'CD4 Th17', 'CD4 Th22', 'CD4 Temra', 'CD4 exhausted',
                                               'CD8 Naive', 'CD8 Tcm CCR4-', 'CD8 Tcm CCR4+', 'CD8 Trm',
                                               'CD8 Tem GZMK+', 'CD8 Tem GZMB+', 'CD8 Temra', 'CD8 NKT-like',
                                               'gdT2', 'MAIT', 'Proliferative T/NK'))
my_palette <- c('#d60000',  '#ff7266', '#edb8b8', '#e36414',
                '#ff9f1c', '#ffbf69', '#f9cb9cff', '#fdf490', 
                '#a3b18a', '#588157', '#3a5a40', '#bce784',
                '#55a630', '#007f5f', '#a5be00', '#74d3ae',
                '#0774d8', '#00acc6','#9651c4')
names(my_palette) <- levels(clon_7449_freq$cell_type)

# 1) plot using stacked barplot
ggplot(clon_7449_freq) +
  geom_bar(aes(x = clonotype_id, y = clone_size, fill = cell_type),
           stat = "identity", color='black') +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 10)) +
  xlab("Clonotype ID") +
  ylab("# clones") +
  scale_fill_manual(name = "Cell type", values = my_palette) 

ggsave(paste0('output/figures/TCR/Freq_barplot/cd4_cd8_clon/cdr3nt_v_j.png'),
       width = 4, height = 3.5,  dpi = 300)
ggsave(paste0('output/figures/TCR/Freq_barplot/cd4_cd8_clon/cdr3nt_v_j.pdf'),
       width = 4, height = 3.5,  dpi = 300)

FeaturePlot(expData, c('EBNA1'))

