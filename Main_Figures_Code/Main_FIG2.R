library("ggsci")
library(ggpubr)
library(tidyverse)

My_Theme2 = theme(strip.background=element_rect(colour="white",fill="white"),legend.text = element_text(size=15), legend.title = element_blank(),legend.position="right",plot.title = element_text(hjust = 0.5, colour = "black", size = 15),
                  axis.title.x = element_text(colour = "black", size = 15),
                  axis.text.x = element_text(colour = "black",size = 15),
                  axis.title.y = element_text(colour = "black",size = 15),axis.text.y = element_text(colour = "black",size = 15),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))

## Fig 2B
primary_cell <- read.delim("https://figshare.com/ndownloader/files/50054727")
primary_cell <- ggplot(primary_cell, aes(x=Primary_Cells, y=Junction_Read)) + geom_bar(stat="identity", fill="#79AF97FF", width=0.5)+ My_Theme2+xlab("Cell Types Expressing Junction Reads") +ylab("CPM") + coord_flip()
