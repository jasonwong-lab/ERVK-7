library("ggsci")
library(ggpubr)
library(tidyverse)

## Figure 2-B - IFN g
IFNg_GSE186169 <- read.delim("https://figshare.com/ndownloader/files/47070703")
comparison <- c("Control","IFNg")
ifng <- ggplot(IFNg_GSE186169, aes(x=Treatment, y=log2TPM, fill=Treatment)) + geom_quasirandom()+stat_compare_means(comparisons = comparison,label = "p.signif",label.y = 2)+ geom_boxplot(width=0.2) + facet_wrap(~HML2_1q22, scales='free')+theme(legend.text = element_text(size=10), legend.title = element_blank(),legend.position="bottom",plot.title = element_text(hjust = 0.5, colour = "black", size = 10),
                                                                                                                                                                                                                                                    axis.title.x = element_text(colour = "black", size = 10), axis.ticks.x =  element_blank(),
                                                                                                                                                                                                                                                    axis.text.x =  element_blank(),
                                                                                                                                                                                                                                                    axis.title.y = element_text(colour = "black",size = 10),axis.text.y = element_text(colour = "black",size = 10),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black")) 
ifng <- ifng + scale_fill_manual(values = c("#374E55FF","#6A6599FF"))+ylab("log2 (TPM)")+theme(strip.background=element_rect(colour="black",fill="white"))+ scale_y_continuous(limits=c(0,4.3))

#### Figure 2-C - TNF a
TNFa_GSE34329 <- read.delim("https://figshare.com/ndownloader/files/47070706")
TNFa_GSE34329$Treatment <- gsub("_","\n", TNFa_GSE34329$Treatment, fixed=TRUE)

tnfa <- ggplot(TNFa_GSE34329, aes(x=Treatment, y=log2TPM, fill=Treatment)) + geom_quasirandom()+stat_compare_means(label = "p.signif",ref.group = ".all.",  method="t.test")+ geom_boxplot(width=0.4) + facet_wrap(~HML2_1q22, scales='free')+theme(legend.text = element_text(size=10), legend.title = element_blank(),legend.position="bottom",plot.title = element_text(hjust = 0.5, colour = "black", size = 10),
                                                                                                                                                                                                                                                    axis.title.x = element_text(colour = "black", size = 10), axis.ticks.x =  element_blank(),
                                                                                                                                                                                                                                                    axis.text.x =  element_blank(),
                                                                                                                                                                                                                                                    axis.title.y = element_text(colour = "black",size = 10),axis.text.y = element_text(colour = "black",size = 10),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black")) 
tnfa <- tnfa + scale_fill_jama()+ylab("log2 (TPM)")+theme(strip.background=element_rect(colour="black",fill="white"))+ scale_y_continuous(limits=c(0,3.3))
