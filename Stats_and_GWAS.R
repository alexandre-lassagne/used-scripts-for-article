### Model on phenotypic data ###

tab_effet <- read.csv("/media/alex/Disque/article2/version_finales/data/tab_ml.csv",  sep = ";", dec = ".", head = TRUE)

LMM = lm(log1p(MicroconidiamL) ~ Box + Batch + ID , data = tab_effect)
hist(LMM$residuals)
x11();par(mfrow=c(2,2));
shapiro.test(LM$residuals)
anova(LMM) 
LM.emm <- as.data.frame(lsmeans(LMM, ~ ID))

write.table(LM.emm, file = "/media/alexandre/Disque/Genomes/GWAS/GWAS_male/pheno.corrected.csv", sep = ";", dec = ".",row.names = F )


### MLMM GWAS ###

myG <- read.table("GWAS_male.hapmap", head = F)
myYa  <- read.csv("pheno.corrected.csv", sep = ";", dec = ".", head = TRUE)
myYa$Mat <- as.factor(myYa$Mat)
K <- read.table("kinship.txt", skip = 3, header = F)
myG$V3 = as.numeric(myG$V3)
setwd("/media/alexandre/Disque/Genomes/GWAS/GWAS_male2/")



myGAPITlsmean <- GAPIT(
  Y=myYa[,c(1,2)],
  G=myG,
  KI=K,
  PCA.total=3,
  model = c("MLMM"))


### Local Score on MMLMM gapit ###

setwd("/media/alexandre/Disque/Genomes/GWAS/GWAS_male2/MLMM_gapit/")

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
