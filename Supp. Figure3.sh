### PopLDdecay installation

 git clone https://github.com/hewm2008/PopLDdecay.git 
        cd PopLDdecay; chmod 755 configure; ./configure;
        make;
        mv PopLDdecay  bin/;    #     [rm *.o]

### VCF diploidisation (mandatory for PopLDdecay)

less SNPCalling.vcf | sed "s|\t0:|\t0\/0:|g" | sed 's|\t1:|\t1\/1:|g' | sed 's|\t2:|\t2\/2:|g' | sed 's|\t3:|\t3\/3:|g' | sed 's|\t.:|\t.\/.:|g' | cat > SNPCalling.vcf_diplo.vcf

### production of the Supp. Figure 3

./bin/PopLDdecay -InVCF SNPCalling.vcf_diplo.vcf -MAF 0cd .05 -OutStat Stat_PopLDdecay
