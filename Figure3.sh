
### LDBlockShow installation

git clone https://github.com/hewm2008/LDBlockShow.git
        cd LDBlockShow ; chmod 755 configure  ;  ./configure;
        make;
        mv LDBlockShow  bin/;    #     [rm *.o]

### VCF diploidisation (mandatory for LDBlockShow)

less GWAS_femmat1_filt.vcf | sed "s|\t0:|\t0\/0:|g" | sed 's|\t1:|\t1\/1:|g' | sed 's|\t2:|\t2\/2:|g' | sed 's|\t3:|\t3\/3:|g' | sed 's|\t.:|\t.\/.:|g' | cat > GWAS_femmat1_filt_diplo.vcf




LDBlockShow -InVCF GWAS_male_diplo.vcf -OutPut /media/alexandre/Disque/Genomes/GWAS/GWAS_male/Figures/localpeaksca4_2 
-Cutline 6.75 -Region Scaffold_4:3617045-4128509 -SeleVar 2 -InGWAS /media/alexandre/Disque/Genomes/GWAS/GWAS_male/tab_gwas_corrected/LDB_input_MLMloc4.txt 
-InGFF /media/alexandre/Disque/Genomes/GWAS/GWAS_male/GUY11_PacBio_merge.gff3
