library(FactoMineR)
library(factoextra)

### Model on phenotypic data ###

tab_effect <- read.csv("tab_ml.csv",  sep = ";", dec = ".", head = TRUE)

LMM = lm(log1p(MicroconidiamL) ~ Box + Batch + ID , data = tab_effect)
hist(LMM$residuals)
x11();par(mfrow=c(2,2));
shapiro.test(LM$residuals)
anova(LMM) 
LM.emm <- as.data.frame(lsmeans(LMM, ~ ID))

write.table(LM.emm, file = "pheno.corrected.csv", sep = ";", dec = ".",row.names = F )


### MLMM GWAS ###

myG <- read.table("GWAS_male.hapmap", head = F)
myYa  <- read.csv("pheno.corrected.csv", sep = ";", dec = ".", head = TRUE)
myYa$Mat <- as.factor(myYa$Mat)
K <- read.table("kinship.txt", skip = 3, header = F)
myG$V3 = as.numeric(myG$V3)

myGAPITlsmean <- GAPIT(
  Y=myYa[,c(1,2)],
  G=myG,
  KI=K,
  PCA.total=3,
  model = c("MLMM"))


### Local Score on MMLMM gapit ###


GWAS = read.csv('GAPIT.Association.GWAS_Results.MLMM.lsmean.csv', sep = ",", dec = ".")
GWAS= GWAS[-c(14801,14802),]


mydata = GWAS[,c(2,3,4)]
head(mydata)
dim(mydata)
colnames(mydata)=c('chr', 'pos','pval')
head(mydata)
dim(mydata)
mydata$chr = as.numeric(mydata$chr)
tapply(mydata$chr,mydata$chr, length)
mydata = as.data.table(mydata)
setkey(mydata, chr)

hist(mydata$pval)


Nchr=length(mydata[,unique(chr)])
chrInfo=mydata[,.(L=.N,cor=autocor(pval)),chr]
chrInfo
setkey(chrInfo,chr)
data.table(chr=mydata[,unique(chr),], S=cumsum(c(0,chrInfo$L[-Nchr]))) %>% 
  setkey(.,chr) %>% mydata[.,posT:=pos+S]


qplot(pval, data=mydata, geom='histogram', binwidth=0.1, main='P-values histogram')

hist(mydata$pval)
z=hist(mydata$pval)
chisq.test(z$counts)

qplot(-log10(pval), data=mydata, geom='histogram', binwidth=0.1, main='-log10(p-values) histogram')
mean(-log10(mydata$pval))
max(-log10(mydata$pval))
summary(-log10(mydata$pval))
quantile(-log10(mydata$pval),0.98)

dim(mydata[-log10(mydata$pval)>1.7,])


xi=2

mydata[,score:= -log10(pval)-xi]

# The score mean must be negative
mean(mydata$score)
max(mydata$score)
mydata[,lindley:=lindley(score),chr]
max(mydata$lindley)
min(mydata$lindley)


##

chrInfo[,th:=thresUnif(L, cor, 2, alpha = 0.01),chr]
chrInfo
mydata=mydata[chrInfo]
mydata$test = 10^-mydata$lindley
##


### Multiple Correspondence Analysis on significant SNP from Local Score Result ###

yolo = read.table(file = "TableSNP.txt" , header = T, row.names = 1)

res.mca <- MCA(yolo, 
               ncp = 4, # nbr de composante 
               graph=FALSE)

fviz_screeplot (res.mca, addlabels = TRUE, ylim = c (0, 70))


res.hcpc <- HCPC(res.mca, graph = FALSE, max = 3)

fviz_dend(res.hcpc, show_labels = FALSE)

# Individus
map_clust =fviz_cluster(res.hcpc, geom = "point", main = "Factor map")

eig.val <- res.mca$eig
barplot(eig.val[, 2], 
        names.arg = 1:nrow(eig.val), 
        main = "Variances Explained by Dimensions (%)",
        xlab = "Principal Dimensions",
        ylab = "Percentage of variances",
        col ="steelblue")

write.table(res.hcpc$data.clust, file = "Cluster_marker_regscaff10.txt", quote = F, row.names = T, col.names = T)

### Clustering of strains based on Multi-locus genotypes ###

clustbis = read.table(file = "/Cluster_marker_regscaff10.txt" , header = T, row.names = 1)
pheno_corre = read.csv('pheno.corrected.csv', sep = ";", dec = ".")
pheno_corre$clust = clustbis$clust
pheno_corre$clust = as.factor(pheno_corre$clust)
lmclust = lm(lsmean ~ clust, data = pheno_corre)
hist(lmclust$residuals)
x11();par(mfrow=c(2,2));plot(lmclust)
shapiro.test(lmclust$residuals)
anova(lmclust) 

model= aov(lsmean ~ clust, data = pheno_corre)
TukeyHSD(model, conf.level=.95)
