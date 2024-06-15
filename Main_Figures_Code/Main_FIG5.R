My_Theme = theme(legend.text = element_text(size=10), legend.title = element_blank(),legend.position="bottom",plot.title = element_text(hjust = 0.5, colour = "black", size = 10),
                 axis.title.x = element_text(colour = "black", size = 10),
                 axis.text.x = element_text(colour = "black",size = 10),
                 axis.title.y = element_text(colour = "black",size = 10),axis.text.y = element_text(colour = "black",size = 10),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))

## Load seurat object
## seurat object can be downloaded from below link 
## https://drive.google.com/file/d/1V66VdLsu3LkB63tGK77TFg308nW867Zy/view?usp=drive_link
## czbiohub/scell_lung_adenocarcinoma explains clearly about how they subtypes the cell annotations
## Reference: Maynard A, McCoach CE, Rotow JK, Harris L, Haderk F, Kerr DL, Yu EA, Schenk EL, Tan W, Zee A et al. 2020. Therapy-Induced Evolution of Human Lung Cancer Revealed by Single-Cell RNA Sequencing. Cell 182: 1232-1251 e1222.
load('IM06_ALL_cells_Immune_AND_nonImmune_annotations-004.RData')

## Brie counts for HML2_1q22.long and HML2_1q22.short
## X0: ambiguous counts - ignored
## X1: Specific to HML2_1q22.long
## X2: Specific to HML2_1q22.short
## X3: Can be applied for both HML2_1q22.long and HML2_1q22.short meaning HML2_1q22 region
summary_table <- read.delim("https://figshare.com/ndownloader/files/47071363")

tiss_subset@meta.data$long_exp <- unlist(summary_table[match(rownames(tiss_subset@meta.data),summary_table$cell_barcode),3])
tiss_subset@meta.data$short_exp <- unlist(summary_table[match(rownames(tiss_subset@meta.data),summary_table$cell_barcode),4])
tiss_subset@meta.data$hml2_exp <- unlist(summary_table[match(rownames(tiss_subset@meta.data),summary_table$cell_barcode),5])
tiss_subset@meta.data$short_exp[is.na(tiss_subset@meta.data$short_exp)] <- 0
tiss_subset@meta.data$hml2_exp[is.na(tiss_subset@meta.data$hml2_exp)] <- 0
tiss_subset@meta.data$long_exp[is.na(tiss_subset@meta.data$long_exp)] <- 0

## Figure 5A : UMAPs
tiss_subset@meta.data$own_annot <- tiss_subset@meta.data$general_annotation1
tiss_subset@meta.data$own_annot[tiss_subset@meta.data$Final_annotation=="tumor"] <- "tumor"

full <- DimPlot(object = tiss_subset,group.by="own_annot")+ggtitle("") + xlab("tSNE1") + ylab("tSNE2")+ My_Theme+scale_color_manual(values=c(endothelial="#374E55FF",epithelial="#DF8F44FF",fibroblast="#00A1D5FF",hepatocyte="#80796BFF",immune="#79AF97FF",melanocytes="#6A6599FF",tumor="#B24745FF"))
hml2 <- DimPlot(object = subset(x = tiss_subset, subset = hml2_exp >0),group.by="own_annot")+My_Theme+ggtitle("") + xlab("tSNE1") + ylab("tSNE2")+scale_color_manual(values=c(endothelial="#374E55FF",epithelial="#DF8F44FF",fibroblast="#00A1D5FF",hepatocyte="#80796BFF",immune="#79AF97FF",melanocytes="#6A6599FF",tumor="#B24745FF"))
long <- DimPlot(object = subset(x = tiss_subset, subset = long_exp >0),group.by="own_annot")+My_Theme+ggtitle("") + xlab("tSNE1") + ylab("tSNE2")+scale_color_manual(values=c(endothelial="#374E55FF",epithelial="#DF8F44FF",fibroblast="#00A1D5FF",hepatocyte="#80796BFF",immune="#79AF97FF",melanocytes="#6A6599FF",tumor="#B24745FF"))
short <- DimPlot(object = subset(x = tiss_subset, subset = short_exp >0),group.by="own_annot")+My_Theme+ggtitle("") + xlab("tSNE1") + ylab("tSNE2")+scale_color_manual(values=c(endothelial="#374E55FF",epithelial="#DF8F44FF",fibroblast="#00A1D5FF",hepatocyte="#80796BFF",immune="#79AF97FF",melanocytes="#6A6599FF",tumor="#B24745FF"))

## Figure 5B: cell type proprotions that are expressing HML2_1q22.long in LUAD primary patients cell types
tiss_subset_cell <- subset(x=tiss_subset, subset = histolgy == "Adenocarcinoma" & primary_or_metastaic=="Primary")
cell_prop <- data.frame(subset(x = tiss_subset, subset = long_exp >0 | hml2_exp > 0| short_exp > 0)@meta.data)
cell_prop <- cell_prop[,c(61:64)]
cell_prop_long <- data.frame(cell_prop[cell_prop$long_exp>0,] %>% group_by(own_annot) %>% summarise(Freq=n()))
cell_prop_long$Freq <- cell_prop_long$Freq*100/sum(cell_prop_long$Freq)
cell_prop_long$type <- "HML2_1q22.long"
cell_prop_short <- data.frame(cell_prop[cell_prop$short_exp>0,] %>% group_by(own_annot) %>% summarise(Freq=n()))
cell_prop_short$Freq <- cell_prop_short$Freq*100/sum(cell_prop_short$Freq)
cell_prop_short$type <- "HML2_1q22.short"
cell_prop_hml2 <- data.frame(cell_prop[cell_prop$hml2_exp>0,] %>% group_by(own_annot) %>% summarise(Freq=n()))
cell_prop_hml2$Freq <- cell_prop_hml2$Freq*100/sum(cell_prop_hml2$Freq)
cell_prop_hml2$type <- "HML2_1q22"
combined <- rbind(cell_prop_long, cell_prop_short)
combined <- rbind(combined, cell_prop_hml2)

cell_prop <- ggplot(combined, aes(x=type, y=Freq,fill=own_annot))+xlab("Expressing")+coord_flip() + geom_bar(stat="identity", width=0.6)+ ylab("Prop (%)")+ My_Theme+scale_fill_manual(values=c(endothelial="#374E55FF",epithelial="#DF8F44FF",fibroblast="#00A1D5FF",hepatocyte="#80796BFF",immune="#79AF97FF",melanocytes="#6A6599FF",tumor="#B24745FF"))

## Figure 5C: Proportion of cells expressing HML2_1q22.long for each cell types
cell_prop <- data.frame(tiss_subset@meta.data)
cell_prop$long_exp_mark <- "no_exp"
cell_prop$long_exp_mark[cell_prop$long_exp>0] <- "exp"
cell_prop_df <- cell_prop %>% group_by(Final_annotation, long_exp_mark) %>% summarize(Freq=n())
cell_prop_df <- data.frame(cell_prop_df[!str_detect(cell_prop_df$Final_annotation, "unknown") & !str_detect(cell_prop_df$Final_annotation, "Unknown") ,])
cell_prop_df_exp <- cell_prop_df[cell_prop_df$long_exp_mark=="no_exp",]
temp <- cell_prop_df[cell_prop_df$long_exp_mark!="no_exp",]
cell_prop_df_exp$exp_freq <- unlist(temp[match(cell_prop_df_exp$Final_annotation, temp$Final_annotation),3])
cell_prop_df_exp[is.na(cell_prop_df_exp)] <- 0
cell_prop_df_exp$prop  <- 0
for(i in 1:nrow(cell_prop_df_exp)){
  cell_prop_df_exp$prop[i] <- cell_prop_df_exp$exp_freq[i]*100/sum(cell_prop_df_exp$Freq[i], cell_prop_df_exp$exp_freq[i])
}
cell_prop_df_exp <- cell_prop_df_exp[order(cell_prop_df_exp$prop, decreasing=TRUE),]
cell_prop_df_exp$sum <- rowSums(cell_prop_df_exp[,3:4])
cell_prop_df_exp$Final_annotation <- factor(cell_prop_df_exp$Final_annotation, levels = cell_prop_df_exp$Final_annotation)
cell_prop_df_exp$own_annot <- unlist(cell_prop[match(cell_prop_df_exp$Final_annotation, cell_prop$Final_annotation),"own_annot"])
cell_prop <- ggplot(cell_prop_df_exp[cell_prop_df_exp$sum>600 | cell_prop_df_exp$Final_annotation=="neuro.e cell",], aes(x=Final_annotation, y=prop, fill=own_annot)) + geom_bar(stat="identity") + My_Theme + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_fill_manual(values=c(endothelial="#374E55FF",epithelial="#DF8F44FF",fibroblast="#00A1D5FF",hepatocyte="#80796BFF",immune="#79AF97FF",melanocytes="#6A6599FF",tumor="#B24745FF"))
cell_prop <- cell_prop + xlab("") + ylab("HML2_1q22.long\nexpressing Prop (%)")+coord_flip()

## Figure 5D: GSEA analysis between tumor cells expressing high (read > 10) HML2_1q22.long vs low (read <= 10) or none expressing HML2_1q22.short
scrna_gsea_raw <- read.delim("https://figshare.com/ndownloader/files/47071360")
scrna_gsea_raw <- scrna_gsea_raw[order(abs(scrna_gsea_raw$normalizedEnrichmentScore),decreasing=TRUE),]
scrna_gsea_nec <- scrna_gsea_raw[1:10,]
scrna_gsea_nec <- scrna_gsea_nec[order((scrna_gsea_nec$normalizedEnrichmentScore)),]
scrna_gsea_nec$geneSet <- factor(scrna_gsea_nec$geneSet, levels=scrna_gsea_nec$geneSet)
scrna_gsea_nec$label <- "Neg"
scrna_gsea_nec$label[scrna_gsea_nec$enrichmentScore > 0] <- "Pos"
scrna_gsea <- ggplot(scrna_gsea_nec, aes(x=geneSet, y=normalizedEnrichmentScore, fill=label)) + geom_bar(stat="identity")+coord_flip()+scale_fill_manual(values=c("#91D1C2B2","#7E6148B2"))+My_Theme + ylab("Enrichment Score") + xlab("") 
