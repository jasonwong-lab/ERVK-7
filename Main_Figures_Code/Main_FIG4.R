My_Theme = theme(legend.text = element_text(size=10), legend.title = element_blank(),legend.position="bottom",plot.title = element_text(hjust = 0.5, colour = "black", size = 10),
                 axis.title.x = element_text(colour = "black", size = 10),
                 axis.text.x = element_text(colour = "black",size = 10),
                 axis.title.y = element_text(colour = "black",size = 10),axis.text.y = element_text(colour = "black",size = 10),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))


## Differentially expressed genes between normal vs tumor
resOrdered_nvst_deseq <- read.delim("https://figshare.com/ndownloader/files/47070856")
## Figure 4-B: H3K4me1 and H3K4me3 related genes
TCGA_LUNG_histone <- read.delim("https://figshare.com/ndownloader/files/47070853")
resOrdered_nvst_deseq$up <- "Up in\nNormal"
resOrdered_nvst_deseq$up[resOrdered_nvst_deseq$log2FoldChange>0] <- "Up in\nTumor"
resOrdered_nvst_deseq$type <- unlist(TCGA_LUNG_histone[match(rownames(resOrdered_nvst_deseq), TCGA_LUNG_histone$gene),"type"])
resOrdered_nvst_deseq$histone <- unlist(TCGA_LUNG_histone[match(rownames(resOrdered_nvst_deseq), TCGA_LUNG_histone$gene),"histone"])
resOrdered_nvst_deseq$label <- paste(resOrdered_nvst_deseq$type, resOrdered_nvst_deseq$histone, sep="\n")
histone <- ggplot(data=resOrdered_nvst_deseq[rownames(resOrdered_nvst_deseq) %in% TCGA_LUNG_histone$gene,], aes(x=log2FoldChange, y=-log10(padj), color=label,label=symbol)) +
  geom_point() + theme_minimal() +geom_text_repel(color="black") +My_Theme+ scale_color_jama()


## Figure 4-C: HSAEC and A549 expression
HSAEC_A549 <- read.delim("https://figshare.com/ndownloader/files/47070859")
HSAEC_A549$batch <- HSAEC_A549$treatment
HSAEC_A549$batch[HSAEC_A549$srr %in% c("SRR7687755","SRR7687756","SRR7687757")] <- "A549\nBatch1"
HSAEC_A549$batch[HSAEC_A549$srr %in% c("SRR8245170","SRR8245171","SRR8245172")] <- "A549\nBatch2"
hsaec_a549_long <- ggplot(HSAEC_A549[HSAEC_A549$name=="HML2_1q22.long",], aes(x=treatment, y=count, fill=treatment)) + geom_bar(position = position_dodge2(preserve = "single"), stat = "identity") + xlab("Cell line") + ylab("log2 (HML2_1q22.long TPM)") + My_Theme + scale_fill_manual(values=c("#6A6599FF","#80796BFF"))

## Figure 4-D: normal and tumor organoid expression
organoid <- read.delim("https://figshare.com/ndownloader/files/47070847")
organoid$label <- gsub("_","\n", organoid$label, fixed=TRUE)
organoid_graph <- ggplot(organoid[organoid$name=='HML2_1q22.long',], aes(x=label, y=count, fill=label)) + geom_boxplot(width=0.2)+xlab(" ") + scale_fill_manual(values=c("#B24745FF","#00A1D5FF"))+ylab("log2 (HML2_1q22.long TPM)")
organoid_graph <- organoid_graph+My_Theme