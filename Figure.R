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

### Figure 1 ###

BARstrain<- ggplot(data=myYa, aes(x=reorder(Strain,log(meanBOX)), y=log(meanBOX), fill=ordre)) + geom_bar(stat="identity",position_dodge(), width=0.5) +
  ylab("log(1p of number of microconidia per ml)") + xlab("Strains")+ theme(plot.title = element_text(hjust = 1)) +  theme(plot.title = element_blank(),
        
        axis.text.x = element_text(color = "black", size = 13, face = "bold",angle = 70, hjust = ),
        axis.title.x =element_text(color = "black", size = 13, face = "bold",angle = 0, hjust = ),
        axis.text.y = element_text(color = "black", size = 15, face = "bold",angle = 0, hjust = ),
        axis.title.y = element_text(color = "black", size = 15, face = "bold",angle = 90, hjust = )
  ) 

### Figure 2 ###

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

### Figure 4 ###

clustered_plot = ggplot(pheno_corre, aes(x=clust, y=lsmean, fill=clust)) +
  geom_jitter(position=position_jitter(0.05)) +
  geom_boxplot(alpha=0.4) +
  stat_summary(fun=mean, geom="point", shape=20, size=3, color="red", fill="red") +
  theme_classic()+
  theme(legend.position="none") +
  ylab("log (1 + number of microconidia per ml)")+
  xlab("Haplotype group")+
  scale_fill_manual(values=c("red", "green","#3388DE"))+
  geom_text(aes(x = 1, y = max(pheno_corre$lsmean) + 1.2, label = "a"), size = 4) +
  geom_text(aes(x = 2, y = max(pheno_corre$lsmean) - 1.4, label = "a b"), size = 4) +
  geom_text(aes(x = 3, y = max(pheno_corre$lsmean) + 0.4, label = "b"), size = 4) +
  theme(plot.title = element_blank(),
                                                              
                                                              axis.text.x = element_text(color = "black", size = 12, face = "bold",angle = 0, hjust = ),
                                                              axis.title.x =element_text(color = "black", size = 15, face = "bold",angle = 0, hjust = , vjust = -0.5),
                                                              axis.text.y = element_text(color = "black", size = 12, face = "bold",angle = 0, hjust = ),
                                                              axis.title.y = element_text(color = "black", size = 15, face = "bold",angle = 90, hjust = , vjust = 2)) +theme(plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), "cm"))


map_clust=fviz_cluster(res.hcpc, geom = "point", main = "", show.clust.cent = F) + theme_classic() +
  theme(plot.title = element_blank(),
        
        axis.text.x = element_text(color = "black", size = 12, face = "bold",angle = 0, hjust = ),
        axis.title.x =element_text(color = "black", size = 15, face = "bold",angle = 0, hjust = , vjust = -0.5),
        axis.text.y = element_text(color = "black", size = 12, face = "bold",angle = 0, hjust = ),
        axis.title.y = element_text(color = "black", size = 15, face = "bold",angle = 90, hjust = , vjust = 2),
        legend.text = element_text(size = 15),
        legend.title = element_text(face = "bold", labs(fill="test")),
        #legend.title = element_text(size = 15,face = "bold", hjust = 0.5 ),
        legend.background = element_rect(color = "grey"),
        legend.margin = margin(0.2,0.3,0.2, 0.2, "cm"),
        legend.box.margin = margin(0, 1, 0, 0, "cm")
  ) + 
  theme(plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), "cm"))

clu = plot_grid(map_clust,clustered_plot, ncol = 2, align = "v", rel_widths = c(4, 3),labels = "AUTO")


### Supp. Figure 5 ###

BOXstrain<- ggplot(data=tab_effect, aes(x=Strain, y=log1p(MicroconidiamL), fill = rep)) + geom_boxplot() + facet_grid(.~Batch,  scales='free')+
  ylab("log (1 + number of microconidia per ml") + xlab("Strains")+ theme(plot.title = element_text(hjust = 1)) + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.001)+ geom_jitter(shape=16, position=position_jitter(0.2))+ 
  theme(plot.title = element_blank(),
        
        axis.text.x = element_text(color = "black", size = 7, face = "bold",angle = 90, hjust = ),
        axis.title.x =element_text(color = "black", size = 13, face = "bold",angle = 0, hjust = ),
        axis.text.y = element_text(color = "black", size = 15, face = "bold",angle = 0, hjust = ),
        axis.title.y = element_text(color = "black", size = 15, face = "bold",angle = 90, hjust = )
  ) 


### Export ###


ggsave(BARstrain, file="figure1.svg", width=12, height=9, dpi = 600)
ggsave(figfin, file="figure2.svg", width=12, height=9, dpi = 600)
ggsave(clu, file="figure4.svg", width=12, height=9, dpi = 600)
ggsave(BOXstrain, file="suppfigure5.svg", width=12, height=9, dpi = 600)








