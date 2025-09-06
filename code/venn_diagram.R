library(tidyverse)
library(hrbrthemes)
library(tm)
library(proustr)
library(VennDiagram)

path <- "/storage1/fs1/martyomov/Active/collaborations/romans/Greg_wu/CLOCK/output/"

# upload & merge data
degs_w_p <- read.csv(paste0(path, 'output/DEA/PBMC_W5_10/combined.csv'))
degs_w_p$tissue <- 'PBMC'
degs_w_c <- read.csv(paste0(path, 'output/DEA/CSF_W5_10/combined.csv'))
degs_w_c$tissue <- 'CSF'
all_degs <- bind_rows(degs_w_p, degs_w_c)

# subset DEGs from PBs (CSF & PBMC)
all_degs_pb <- all_degs %>%
  filter(cell_type %in% 'Plasmablast' & p_val_adj < 0.05 & avg_log2FC > 0) %>% 
  rename(gene = V1) %>%
  select(gene, cell_type, tissue) 

# plot
venn.diagram(
  x = list(all_degs_pb %>% filter(tissue == "PBMC") %>% select(gene) %>% unlist() , 
           all_degs_pb %>% filter(tissue == "CSF") %>% select(gene) %>% unlist()),
  category.names = c("PBMC" , "CSF" ),
  filename = paste0(path, 'figures/pdf/venn.png'),
  output = TRUE , imagetype="png" , height = 500 , width = 480 , resolution = 300,
  compression = "lzw", lwd = 1,                              # thickness of the borders
  col=c("#c1121f", '#21908dff'),
  fill = c(alpha("#c1121f", 0.7), alpha('#21908dff', 0.7)),
  cex = 0.5, fontfamily = "sans",    # cex is a font size
  cat.cex = 0.5, cat.default.pos = "outer",
  cat.pos = c(0, 0), cat.dist = c(0.03, -0.425),
  cat.fontfamily = "sans", 
  cat.col = c("#c1121f", '#21908dff'),   # labels color
  rotation.degree = 270,      # rotate circules
  scaled = F,                 # disable scaling --> the same size of circules 
  disable.logging = T)
