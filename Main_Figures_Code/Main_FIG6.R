My_Theme = theme(legend.text = element_text(size=10), legend.title = element_blank(),legend.position="right",plot.title = element_text(hjust = 0.5, colour = "black", size = 10),
                 axis.title.x = element_text(colour = "black", size = 10),
                 axis.text.x = element_text(colour = "black",size = 10),
                 axis.title.y = element_text(colour = "black",size = 10),axis.text.y = element_text(colour = "black",size = 10),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))

My_Theme2 = theme(legend.text = element_text(size=10), legend.title = element_blank(),legend.position="bottom",plot.title = element_text(hjust = 0.5, colour = "black", size = 10),
                 axis.title.x = element_text(colour = "black", size = 10),
                 axis.text.x = element_text(colour = "black",size = 10),
                 axis.title.y = element_text(colour = "black",size = 10),axis.text.y = element_text(colour = "black",size = 10),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))

library("ggsci")
library(ggpubr)
library(tidyverse)
library(org.Hs.eg.db)
library(msigdbr)

## Fig 6B
TNFa_GSE34329 <- read.delim("https://figshare.com/ndownloader/files/49913514")
TNFa_GSE34329$Treatment <- gsub("_","\n", TNFa_GSE34329$Treatment, fixed=TRUE)

tnfa <- ggplot(TNFa_GSE34329, aes(x=Treatment, y=log2TPM, fill=Treatment)) + geom_quasirandom()+stat_compare_means(label = "p.signif",ref.group = ".all.",  method="t.test")+ geom_boxplot(width=0.4) + facet_wrap(~HML2_1q22, scales='free')+theme(legend.text = element_text(size=10), legend.title = element_blank(),legend.position="bottom",plot.title = element_text(hjust = 0.5, colour = "black", size = 10),
                                                                                                                                                                                                                                                    axis.title.x = element_text(colour = "black", size = 10), axis.ticks.x =  element_blank(),
                                                                                                                                                                                                                                                    axis.text.x =  element_blank(),
                                                                                                                                                                                                                                                    axis.title.y = element_text(colour = "black",size = 10),axis.text.y = element_text(colour = "black",size = 10),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black")) 
tnfa <- tnfa + scale_fill_jama()+ylab("log2 (TPM)")+theme(strip.background=element_rect(colour="black",fill="white"))+ scale_y_continuous(limits=c(0,3.3))

## Fig 6C
fig6C <- read.delim("/storage2/jwlab/sandy/LUNG_paper/data/fig6C.txt")
tnfa_drug <- ggplot(fig6C, aes(x=type, y=logval, fill=type, group=type))+scale_y_continuous(limits=c(0,5))+ scale_fill_manual(values=c("#ED0000FF","#42B540FF","#925E9FFF"))+theme(strip.background=element_rect(colour="white",fill="white"))+geom_point(shape=21) + My_Theme2 + facet_wrap(~name)+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("Treatment") + ylab("log2 (TPM)") +stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,  geom = "crossbar",  width = 0.3,color ="#DF8F44FF") 

## Fig 6D
resOrdered_hvsl_deseq <- read.delim("https://figshare.com/ndownloader/files/49913508")
tcga_gsea <- c(resOrdered_hvsl_deseq$log2FoldChange)
names(tcga_gsea) <- rownames(resOrdered_hvsl_deseq)
tcga_gsea <- na.omit(tcga_gsea)
tcga_gsea = sort(tcga_gsea, decreasing = TRUE)
tnf_receptor <- read.table("https://figshare.com/ndownloader/files/49913505", quote="\"", comment.char="")
tnf_receptor$ensembl <-mapIds(x = org.Hs.eg.db,
                              keys = tnf_receptor$V1,
                              "ENSEMBL",
                              "SYMBOL",
                              fuzzy = TRUE)

tnf_receptor$V1 <- "TNF_receptor"
gse <- GSEA(tcga_gsea, TERM2GENE = tnf_receptor, pvalueCutoff =5, exponent=0.3)
#NES = 1.759212 p-value = 0.01106049
tcga_gsea <- gseaplot2(gse, geneSetID = 1, title = "TNF Receptor Signaling Pathway",  color="#80796BFF") +theme_classic()

## Fig 6E
TCGA_LUNG_TPM <- read.delim("https://figshare.com/ndownloader/files/47070841")
luad_pat <- read.delim("https://figshare.com/ndownloader/files/47070835", header=FALSE)
luad_pat <- data.frame(t(luad_pat))
TCGA_LUNG_TPM_nec <- TCGA_LUNG_TPM[,substr(colnames(TCGA_LUNG_TPM), 14, 16) %in% c("01A","11A")]
TCGA_LUNG_TPM_luad <- TCGA_LUNG_TPM_nec[,colnames(TCGA_LUNG_TPM_nec) %in% luad_pat$t.luad_pat.]
rela_exp <- TCGA_LUNG_TPM_luad[rownames(TCGA_LUNG_TPM_luad)=="ENSG00000173039",]
long <- TCGA_LUNG_TPM_luad[1,]
combine <- data.frame(t(rbind(long , rela_exp)))
combine$HML2_1q22.long <- log(combine$HML2_1q22.long+1, base=2)
combine$ENSG00000173039 <- log(combine$ENSG00000173039+1, base=2)
combine$long_type <- "High\n(N=255)"
combine$long_type[combine$HML2_1q22.long < median(combine$HML2_1q22.long)] <- "Low\n(N=255)"
combine$type2 <- NA
combine$type2[combine$HML2_1q22.long > quantile(combine$HML2_1q22.long, 0.75)] <- "High (N=128)"
combine$type2[combine$HML2_1q22.long < quantile(combine$HML2_1q22.long, 0.25)] <- "Low (N=128)"

rela_imp <- ggplot(na.omit(combine), aes(x=type2, y=ENSG00000173039)) + geom_quasirandom(shape=21, fill="white", alpha=0.2) + geom_boxplot(width=0.2) + My_Theme+ xlab("ERVK-7.long Expression") + ylab("log2 (RELA TPM)")  
ggsave(file="/storage2/jwlab/sandy/LUNG_paper/figure/supple/rela_imp.pdf", plot=rela_imp,bg = 'white', width = 8, height = 8, units = 'cm', dpi = 600)


## Fig 6F
IFNg_GSE186169 <- read.delim("https://figshare.com/ndownloader/files/49913511")
comparison <- c("Control","IFNg")
ifng <- ggplot(IFNg_GSE186169, aes(x=Treatment, y=log2TPM, fill=Treatment)) + geom_quasirandom()+stat_compare_means(comparisons = comparison,label = "p.signif",label.y = 2)+ geom_boxplot(width=0.2) + facet_wrap(~HML2_1q22, scales='free')+theme(legend.text = element_text(size=10), legend.title = element_blank(),legend.position="bottom",plot.title = element_text(hjust = 0.5, colour = "black", size = 10),
                                                                                                                                                                                                                                                    axis.title.x = element_text(colour = "black", size = 10), axis.ticks.x =  element_blank(),
                                                                                                                                                                                                                                                    axis.text.x =  element_blank(),
                                                                                                                                                                                                                                                    axis.title.y = element_text(colour = "black",size = 10),axis.text.y = element_text(colour = "black",size = 10),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black")) 
ifng <- ifng + scale_fill_manual(values = c("#374E55FF","#6A6599FF"))+ylab("log2 (TPM)")+theme(strip.background=element_rect(colour="black",fill="white"))+ scale_y_continuous(limits=c(0,4.3))

## Fig 6G
res_select1_df <- read.delim("https://figshare.com/ndownloader/files/49913520")
pathwaysH_tnfa2 <- list(pathwaysH$HALLMARK_INTERFERON_GAMMA_RESPONSE)
pathwaysH_tnfa2 <- data.frame(pathwaysH_tnfa2)
pathwaysH_tnfa2$label <- "Gamma"
colnames(pathwaysH_tnfa2)[1] <- "gene_name"
pathwaysH_tnfa2 <- pathwaysH_tnfa2[,c(2,1)]
pathwaysH_tnfa_alpha <- pathwaysH_tnfa
pathwaysH_tnfa <- rbind(pathwaysH_tnfa, pathwaysH_tnfa2)

tcga_gsea_hist <- c(res_select1_df$log2FoldChange)
names(tcga_gsea_hist) <- rownames(res_select1_df)
tcga_gsea_hist <- na.omit(tcga_gsea_hist)
tcga_gsea_hist = sort(tcga_gsea_hist, decreasing = TRUE)
gse <- GSEA(tcga_gsea_hist, TERM2GENE = pathwaysH_tnfa, pvalueCutoff =5)

tcga_gsea_hist <- gseaplot2(gse, geneSetID = 1, title = "INTERFERON_GAMMA_RESPONSE",  color=c("#6A6599FF")) +theme_classic()

## Fig 6H
ervk7 <- read.delim("https://figshare.com/ndownloader/files/49913517")
ervk7_lusc <- ervk7[ervk7$cohort=="LUSC",]
for(i in c(7:ncol(ervk7_lusc))){
  ervk7_lusc[,i] <- as.numeric(ervk7_lusc[,i])
  ervk7_lusc[,i] <- (ervk7_lusc[,i]-mean(ervk7_lusc[,i]))/sd(ervk7_lusc[,i])
}
mean_cat <- nrow(ervk7_lusc[ervk7_lusc$copynum!="0",])
mean_cat <- mean_cat/nrow(ervk7_lusc)
sd_val <- sqrt(mean_cat*(1-mean_cat))
ervk7_lusc$copynum <- as.numeric(ervk7_lusc$copynum)
ervk7_lusc$copynum <- (ervk7_lusc$copynum - mean_cat)/sd_val

mlm1 <- lm(HML2_1q22 ~  copynum + ifngamma+nfkb, data = ervk7_lusc)

summary(mlm1)

ervk7_luad <- ervk7[ervk7$cohort=="LUAD",]
for(i in c(7:ncol(ervk7_luad))){
  ervk7_luad[,i] <- as.numeric(ervk7_luad[,i])
  ervk7_luad[,i] <- (ervk7_luad[,i]-mean(ervk7_luad[,i]))/sd(ervk7_luad[,i])
}
mean_cat <- nrow(ervk7_luad[ervk7_luad$copynum!="0",])
mean_cat <- mean_cat/nrow(ervk7_luad)
sd_val <- sqrt(mean_cat*(1-mean_cat))
ervk7_luad$copynum <- as.numeric(ervk7_luad$copynum)
ervk7_luad$copynum <- (ervk7_luad$copynum - mean_cat)/sd_val

mlm1_luad <- lm(HML2_1q22 ~ copynum+ ifngamma+nfkb, data = ervk7_luad)

summary(mlm1_luad)

mlm1_df <- (summary(mlm1))
mlm1_df_plot <- mlm1_df$coefficients
lb = confint(mlm1)[,1] # lower bound confidence 
ub = confint(mlm1)[,2]
mlm1_df_plot <- mlm1_df_plot[2:nrow(mlm1_df_plot),]
mlm1_df_plot <- data.frame(mlm1_df_plot)
mlm1_df_plot$lower <- c(lb[2:length(lb)])
mlm1_df_plot$upper <- c(ub[2:length(ub)])
mlm1_df_plot$name <- rownames(mlm1_df_plot)

mlm1_luad_df <- (summary(mlm1_luad))
mlm1_luad_df_plot <- mlm1_luad_df$coefficients
lb = confint(mlm1_luad)[,1] # lower bound confidence 
ub = confint(mlm1_luad)[,2]
mlm1_luad_df_plot <- mlm1_luad_df_plot[2:nrow(mlm1_luad_df_plot),]
mlm1_luad_df_plot <- data.frame(mlm1_luad_df_plot)
mlm1_luad_df_plot$lower <- c(lb[2:length(lb)])
mlm1_luad_df_plot$upper <- c(ub[2:length(ub)])
mlm1_luad_df_plot$name <- rownames(mlm1_luad_df_plot)

combined_model <- rbind(mlm1_df_plot, mlm1_luad_df_plot)
combined_model$cohort <- "LUAD"
combined_model$cohort[1:3] <- "LUSC"


library(ggplot2)
forest_plot <- ggplot(data=combined_model, aes(x=name, y=Estimate, ymin=lower, ymax=upper, color=cohort, group=cohort)) +
  geom_pointrange(aes(group=cohort), position=position_dodge(width=0.3)) +
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +My_Theme2 + ylab("Odd Ratio")+ xlab("") + scale_color_manual(values=c("#B24745FF","#00a1d5ff"))

