
### LDBlockShow installation

git clone https://github.com/hewm2008/LDBlockShow.git
        cd LDBlockShow ; chmod 755 configure  ;  ./configure;
        make;
        mv LDBlockShow  bin/;    #     [rm *.o]

### VCF diploidisation (mandatory for LDBlockShow)

less SNPCalling.vcf | sed "s|\t0:|\t0\/0:|g" | sed 's|\t1:|\t1\/1:|g' | sed 's|\t2:|\t2\/2:|g' | sed 's|\t3:|\t3\/3:|g' | sed 's|\t.:|\t.\/.:|g' | cat > SNPCalling.vcf_diplo.vcf




LDBlockShow -InVCF SNPCalling.vcf_diplo.vcf -OutPut Figure3 
-Cutline 4.20 -Region Scaffold_10:917290-968350 -SeleVar 2 -InGWAS /media/alexandre/Disque/Genomes/GWAS/GWAS_male/tab_gwas_corrected/LDB_input_MLMloc4.txt 
-InGFF GUY11_PacBio_merge.gff3
