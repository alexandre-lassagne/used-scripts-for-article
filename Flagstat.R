library(ggplot2)
library(gridExtra)
library(ggsave)

propmapA <- read.csv("propmapped_final.csv", sep = ";", header = T, dec = ".")
depthA <- read.csv("DP_MA_final.csv",  sep = ";", header = T, dec = ".")

### Percentage of mapped reads on the reference genome 

test2 <- ggplot(data=propmapA, aes(x=ID, y=PROPORTION)) + geom_bar(stat="identity", width=0.7)+ theme(legend.position="none") + 
  geom_hline(yintercept=75, linetype="dashed", color = "red") +
  ylab("percentage of read mapped") + xlab("Strains")+ theme(plot.title = element_text(hjust = 1)) + 
  ylim(0, 100) +
  theme(plot.title = element_blank(),
        
        axis.text.x = element_text(color = "black", size = 7, face = "bold",angle = 90, hjust = ),
        axis.title.x =element_text(color = "black", size = 13, face = "bold",angle = 0, hjust = ),
        axis.text.y = element_text(color = "black", size = 15, face = "bold",angle = 0, hjust = ),
        axis.title.y = element_text(color = "black", size = 15, face = "bold",angle = 90, hjust = )
  ) +
  theme(
    plot.title = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y.right =element_blank() )


### Number of mapped reads (in blue) on total reads (in yellow)

test3 <- ggplot(data=propmapA) +
  geom_bar(stat="identity", width=0.7, aes(x=ID, y=(TOT/1000000)), color="yellow", fill="yellow") + 
  geom_bar(stat="identity", width=0.7, aes(x=ID, y=(MAPPED_reads/1000000)), color="yellow", fill="blue")+
  theme(legend.position="none") + 
  ylab("number of reads (million)") + xlab("Strains") + theme(plot.title = element_text(hjust = 1)) + 
  theme(plot.title = element_blank(),
        
        axis.text.x = element_text(color = "black", size = 9, face = "bold",angle = 70, hjust = +1),
        axis.title.x =element_text(color = "black", size = 15, face = "bold",angle = 0, hjust = ),
        axis.text.y = element_text(color = "black", size = 15, face = "bold",angle = 0, hjust = ),
        axis.title.y = element_text(color = "black", size = 15, face = "bold",angle = 90, hjust = )
  ) 


### mean reads depth per isolates

test4 <- ggplot(data=depthA, aes(x=ID, y=mean_read_DP)) + geom_bar(stat="identity", width=0.7)+ theme(legend.position="none") + 
  geom_hline(yintercept=63, linetype="dashed", color = "red") +
  ylab("mean depth") + xlab("Strains") + theme(plot.title = element_text(hjust = 1)) + 
  ylim(0, 280) + 
  theme(plot.title = element_blank(),
        
        axis.text.x = element_text(color = "black", size = 7, face = "bold",angle = 90, hjust = ),
        axis.title.x =element_text(color = "black", size = 13, face = "bold",angle = 0, hjust = ),
        axis.text.y = element_text(color = "black", size = 15, face = "bold",angle = 0, hjust = ),
        axis.title.y = element_text(color = "black", size = 15, face = "bold",angle = 90, hjust = )
  ) +
  theme(
    plot.title = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y.right =element_blank() )

### percentage covered on the reference genome 


test5  <- ggplot(data=depthA, aes(x=ID, y=breath_DP)) + geom_bar(stat="identity", width=0.7)+ theme(legend.position="none") + 
  geom_hline(yintercept=80, linetype="dashed", color = "red") +
  ylab("percentage covering") + xlab("Strains") + theme(plot.title = element_text(hjust = 1)) + 
  ylim(0, 100) +
  theme(plot.title = element_blank(),
        
        axis.text.x = element_text(color = "black", size = 9, face = "bold",angle = 70, hjust = +1),
        axis.title.x =element_text(color = "black", size = 15, face = "bold",angle = 0, hjust = ),
        axis.text.y = element_text(color = "black", size = 15, face = "bold",angle = 0, hjust = ),
        axis.title.y = element_text(color = "black", size = 15, face = "bold",angle = 90, hjust = )
  ) 

bigtest <- grid.arrange(test2,test4,test3,test5,ncol=2, nrow = 2)

ggsave(bigtest, file="Suppfigure1.svg", width=19, height=9, dpi = 600)
