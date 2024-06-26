library(vcfR)
library(adegenet)

pop.data <- read.csv("popdata.csv", sep = ";", header = TRUE)
vcf_oryzaeY <- read.vcfR("\GWAS_male_13sca.vcf")

genlight_oryzaeY <- vcfR2genlight(vcf_oryzaeY)

ploidy(genlight_oryzaeY) <- 1
genlight_oryzaeY

all(colnames(genlight_oryzaeY@pop)[-1] == pop.data$population)
pop (genlight_oryzaeY) <- pop.data$population
genlight_oryzaeY




acp = glPca(genlight_oryzaeY, center = TRUE, scale = FALSE, nf = NULL, loadings = TRUE, 
      alleleAsUnit = FALSE, useC = TRUE, parallel = FALSE,
      n.cores = NULL, returnDotProd=FALSE, matDotProd=NULL)

scatter.glPca(acp, posi.da="bottomright", bg="white", pch=17:22)
