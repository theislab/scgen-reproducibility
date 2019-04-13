library(dplyr)
library("ComplexHeatmap")
#################
load("./BlueYellowColormaps_V1.RData")
pdf(file="../results/Figures/Supplemental Figure 8/SupplFig8c_salmonella_TA_heatmap.pdf")
heat = read.csv("./TA.txt")
dim(heat)
df = heat[,2:11]
dim(df)
d= df %>% as.matrix() %>% t
d = t(apply(d, 1, scale))
font_size = .6
ha_row = rowAnnotation(df = data.frame(genes = c(rep("upregulated", 5), rep("downregulated", 5))),
                       col = list(genes = c("downregulated" =  "orchid", "upregulated" = "yellowgreen")),
                       show_annotation_name = TRUE,annotation_name_side = "top",annotation_name_rot =00,annotation_name_gp = gpar(cex =.8) )
ht = Heatmap(d,row_names_side = "left",
             heatmap_legend_param = list(title = " scaled expression")
             , cluster_rows = FALSE,cluster_columns = FALSE,column_title = "intermediate TA cells",col = yellow2blue)
ht_list =  ht  + ha_row
draw(ht_list)
dev.off()
#####################
pdf(file="../results/Figures/Supplemental Figure 8/SupplFig8d_pbmc_CD4T_heatmap.pdf")
heat = read.csv("./CD4T.txt")
dim(heat)
df = heat[,2:11]
dim(df)
d= df %>% as.matrix() %>% t
d = t(apply(d, 1, scale))
font_size = .6
ht = Heatmap(d,row_names_side = "left",
             heatmap_legend_param = list(title = "scaled expression")
             , cluster_rows = FALSE,cluster_columns = FALSE,column_title = "intermediate CD4-T cells",col = yellow2blue)
ht_list =  ht
draw(ht_list)
dev.off()

