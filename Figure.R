library(ggmanh)
library(ggplot2)
library(lsmeans)
library(tidyverse)
library(data.table)
library(gridExtra)
library(cowplot)
library(FactoMineR)
library(factoextra)
library(ggsave)



BARstrain<- ggplot(data=myYa, aes(x=reorder(Strain,log(meanBOX)), y=log(meanBOX), fill=ordre)) + geom_bar(stat="identity",position_dodge(), width=0.5) +
  ylab("log(1p of number of microconidia per ml)") + xlab("Strains")+ theme(plot.title = element_text(hjust = 1)) +  theme(plot.title = element_blank(),
        
        axis.text.x = element_text(color = "black", size = 13, face = "bold",angle = 70, hjust = ),
        axis.title.x =element_text(color = "black", size = 13, face = "bold",angle = 0, hjust = ),
        axis.text.y = element_text(color = "black", size = 15, face = "bold",angle = 0, hjust = ),
        axis.title.y = element_text(color = "black", size = 15, face = "bold",angle = 90, hjust = )
  ) 

### Figure GWAS

g <- manhattan_plot(x = mydata , pval.colname = "pval", chr.colname = "chr", pos.colname = "pos", plot.title = "GWAS", signif = 0.000004, y.label = "-log10(p_value)") 
plotgwas = g + ylim(0,5.8) +
  theme(
    plot.title = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y.right =element_blank() ) + theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))

p <- manhattan_plot(x = mydata, pval.colname = "test", chr.colname = "chr", pos.colname = "pos", plot.title = "", y.label = "Local Score", cutoff = NULL, signif = 1e-60, x.label = "Scaffolds")

local_score = p +   geom_errorbar(aes(ymax=th, ymin=th), linetype="dashed", position=position_dodge()) + ylim(0,20.5) + 
  theme(
    plot.title = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y.right =element_blank() )

figfin = plot_grid(plotgwas,local_score, ncol = 1, align = "v", rel_heights=c(3, 6))



ggsave(BARstrain, file="figure1.svg", width=12, height=9, dpi = 600)
ggsave(figfin, file="figure2.svg", width=12, height=9, dpi = 600)
