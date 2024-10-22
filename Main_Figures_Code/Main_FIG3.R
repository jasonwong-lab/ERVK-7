My_Theme = theme(legend.text = element_text(size=10), legend.title = element_blank(),legend.position="bottom",plot.title = element_text(hjust = 0.5, colour = "black", size = 10),
                 axis.title.x = element_text(colour = "black", size = 10),
                 axis.text.x = element_text(colour = "black",size = 10),
                 axis.title.y = element_text(colour = "black",size = 10),axis.text.y = element_text(colour = "black",size = 10),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))

library(tidyverse)
library(ggsci)

TCGA_LUNG_TPM <- read.delim("https://figshare.com/ndownloader/files/47070841")
luad_pat <- read.delim("https://figshare.com/ndownloader/files/47070835", header=FALSE)
luad_pat <- data.frame(t(luad_pat))

##Fig 3C
TCGA_LUNG_TPM_nec <- data.frame(melt(cbind(name=rownames(TCGA_LUNG_TPM)[1:4], TCGA_LUNG_TPM[1:4,])))
TCGA_LUNG_TPM_nec$value <- log(as.numeric(TCGA_LUNG_TPM_nec$value)+1, base=2)
TCGA_LUNG_TPM_nec$cohort <- "LUSC"
TCGA_LUNG_TPM_nec$cohort[TCGA_LUNG_TPM_nec$variable %in% luad_pat$t.luad_pat.] <- "LUAD"
TCGA_LUNG_TPM_nec$label <- unlist(sapply(strsplit(as.character(TCGA_LUNG_TPM_nec$variable), ".", fixed=TRUE), function(x) x[4], simplify=FALSE))
TCGA_LUNG_TPM_nec <- TCGA_LUNG_TPM_nec[TCGA_LUNG_TPM_nec$label %in% c("01A","11A"),]
TCGA_LUNG_TPM_nec$label[TCGA_LUNG_TPM_nec$label=="11A"] <- "Normal"
TCGA_LUNG_TPM_nec$label[TCGA_LUNG_TPM_nec$label=="01A"] <- "Tumor"
TCGA_LUNG_TPM_nec$final_label <- paste(TCGA_LUNG_TPM_nec$cohort,TCGA_LUNG_TPM_nec$label, sep="\n")
TCGA_LUNG_TPM_nec$final_label[str_detect(TCGA_LUNG_TPM_nec$final_label, "Normal")] <- "Normal"
my_comparisons <- list(c("Normal", "LUAD\nTumor"), c("Normal", "LUSC\nTumor"), c("LUAD\nTumor", "LUSC\nTumor"))
TCGA_LUNG_TPM_nec$name[TCGA_LUNG_TPM_nec$name == "HML2_1q22.long"] <- "ERVK-7.long"
TCGA_LUNG_TPM_nec$name[TCGA_LUNG_TPM_nec$name == "HML2_1q22"] <- "ERVK-7"
TCGA_LUNG_TPM_nec$name[TCGA_LUNG_TPM_nec$name == "HML2_1q22.short"] <- "ERVK-7.short"

tcga <- ggplot(TCGA_LUNG_TPM_nec[str_detect(TCGA_LUNG_TPM_nec$name, "ERVK-7"),], aes(x=final_label, y=value, fill=final_label))+ geom_quasirandom(shape=21, fill="white", alpha=0.1) + geom_boxplot(width=0.2,outlier.shape = NA)+ facet_wrap(~name)+ scale_fill_manual(values=c("#B24745FF","#00A1D5FF","#DF8F44FF"))
tcga <- tcga + xlab("") + ylab("log2 (TPM)") + stat_compare_means(comparisons = my_comparisons, label = "p.signif",ref.group = ".all.", method="t.test")+theme(legend.text = element_text(size=10), legend.title = element_blank(),legend.position="bottom",plot.title = element_text(hjust = 0.5, colour = "black", size = 10),
                                                                                                                                                               axis.title.x = element_text(colour = "black", size = 10), axis.ticks.x =  element_blank(),
                                                                                                                                                               axis.text.x =  element_blank(),
                                                                                                                                                               axis.title.y = element_text(colour = "black",size = 10),axis.text.y = element_text(colour = "black",size = 10),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black")) 
tcga <- tcga + scale_y_continuous(limits=c(0,9))+theme(strip.background=element_rect(colour="white",fill="white"))

## correlation
TCGA_LUNG_TPM_cor <- data.frame(t(TCGA_LUNG_TPM[1:4,]))
TCGA_LUNG_TPM_cor$final_label <- unlist(TCGA_LUNG_TPM_nec[match(rownames(TCGA_LUNG_TPM_cor), TCGA_LUNG_TPM_nec$variable),"final_label"])
TCGA_LUNG_TPM_cor <- na.omit(TCGA_LUNG_TPM_cor)
for(i in 1:4){
  TCGA_LUNG_TPM_cor[,i] <- log(as.numeric(TCGA_LUNG_TPM_cor[,i])+1, base=2)
}
cor1 <- ggplot(TCGA_LUNG_TPM_cor, aes(x=HML2_1q22.long, HML2_1q22, color=final_label)) + geom_point(fill="black", color="black",alpha=0.1, shape=21)+My_Theme+ scale_color_manual(values=c("#B24745FF","#00A1D5FF","#DF8F44FF"))
cor1 <- cor1 + xlab("log2 (ERVK-7.long TPM)") + ylab("log2 (ERVK-7 TPM)")+ facet_wrap(~final_label, scales="free")
cor1 <- cor1 + stat_cor(method = "spearman", hjust=0.5, vjust=0, label.x = TCGA_LUNG_TPM_cor$HML2_1q22.long, 
                        label.y=TCGA_LUNG_TPM_cor$HML2_1q22)+stat_smooth(method='lm', formula = 'y ~ x',size=0.8)+ scale_x_continuous(limits=c(0,7))+ scale_y_continuous(limits=c(0,8))
cor2 <- ggplot(TCGA_LUNG_TPM_cor, aes(x=HML2_1q22.short, HML2_1q22, color=final_label)) + geom_point(fill="black",color="black" ,alpha=0.1, shape=21)+My_Theme+ scale_color_manual(values=c("#B24745FF","#00A1D5FF","#DF8F44FF"))
cor2 <- cor2 + xlab("log2 (ERVK-7.short TPM)") + ylab("log2 (ERVK-7 TPM)")+ facet_wrap(~final_label, scales="free")
cor2 <- cor2 + stat_cor(method = "spearman", hjust=0.5, vjust=0, label.x = TCGA_LUNG_TPM_cor$HML2_1q22.long, 
                        label.y=TCGA_LUNG_TPM_cor$HML2_1q22)+stat_smooth(method='lm', formula = 'y ~ x',size=0.8)+ scale_x_continuous(limits=c(0,8.5))+ scale_y_continuous(limits=c(0,8))
cor1 <- cor1 + theme(strip.background=element_rect(colour="white",fill="white"))
cor2 <- cor2 + theme(strip.background=element_rect(colour="white",fill="white"))
