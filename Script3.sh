
#!/bin/sh

# Function to display usage message
usage() {
    echo "Usage: $0 
        This script allows to perfom SNPCalling anf SNPFiltering
        It will create 2 outpout: 
        an unfiltered compressed vcf
        a filtered vcf
        
        default values:
        THREADS=6
        MINDP=10 
        NONMISS=0.9
        MAF=0.05
        
        indels are removed but polyallelic SNP are kept 

        
        Flags option:
    [-t threads] 
    -i input_mapping
    -f fastq_path 
    -o output_path
    -m minimum depth
    -n percentage of non missing data
    -M Minimum Allele Frequencies 
    "
    exit 1
}

# Default values
THREADS=6
MINDP=10 
NONMISS=0.9
MAF=0.05


# Parse command-line arguments
while getopts t:i:r:o:m:n:M:h flag
do
    case "${flag}" in
        t) THREADS=${OPTARG};;
        i) INPUT=${OPTARG};;
        r) REF=${OPTARG};;
        o) OUTPUT=${OPTARG};;
        m) MINDP=${OPTARG};;
        n) NONMISS=${OPTARG};;
        M) MAF=${OPTARG};;
        h) usage;;
        *) usage;;
    esac
done

# Check if all required parameters are provided
if [ -z "$REF" ] || [ -z "$OUTPUT" ] || [ -z "$INPUT" ]; then
    echo "Error: Missing required arguments."
    usage
fi

bcftools mpileup --threads 20 --annotate "DP,AD,INFO/AD" -f ${REF} ${INPUT}/*.filt_sort.bam | bcftools call --threads 20 -m --ploidy 1 -f GQ| bgzip > ${OUTPUT}/SNPCalling_unfilt.vcf.gz

vcftools --gzvcf ${OUTPUT}/SNPCalling_unfilt.vcf.gz --remove-indels --minDP ${MINDP} --max-missing ${NONMISS} --maf ${MAF} --recode --out ${OUTPUT}/SNPCalling_minDP${MINDP}_nonmiss${NONMISS}_maf${MAF}
