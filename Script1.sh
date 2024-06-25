#!/bin/sh

# Function to display usage message
usage() {
    echo "Usage: $0 
        This script allows to perfom mapping with bwa mem and samtools
        It will create 2 files in the output_path: 
        Mapping => reads mapped
        NOMapping => reads non-mapped
        
        Flags option:
    -t threads 
    -r reference 
    -f fastq_path 
    -o output_path
    -q mapping quality
    "
    exit 1
}

# Parse command-line arguments
while getopts t:r:f:o:q:h flag
do
    case "${flag}" in
        t) THREADS=${OPTARG};;
        r) REF=${OPTARG};;
        f) PATH_FASTQ=${OPTARG};;
        o) OUTPUT=${OPTARG};;
        q) QUAL=${OPTARG};;
        h) usage;;
        *) usage;;
    esac
done

# Check if all required parameters are provided
if [ -z "$THREADS" ] || [ -z "$REF" ] || [ -z "$PATH_FASTQ" ] || [ -z "$OUTPUT" ] || [ -z "$QUAL" ]; then
    echo "Error: Missing required arguments."
    usage
fi

# Construct output directories
OUTPUT_MAPPING="${OUTPUT}/Mapping/"
OUTPUT_NOMAPPING="${OUTPUT}/NOMapping/"
mkdir -p "$OUTPUT"
mkdir -p "$OUTPUT_MAPPING"
mkdir -p "$OUTPUT_NOMAPPING"

## Indexing the reference genome
bwa index -a bwtsw $REF

for fastq in ${PATH_FASTQ}*R1.fastq.gz; do 


  ID=`basename -s _R1.fastq.gz $fastq`

  bwa mem -t ${THREADS} $REF ${PATH_FASTQ}${ID}_R1.fastq.gz ${PATH_FASTQ}${ID}_R2.fastq.gz -R "@RG\tID:${ID}\tSM:${ID}\tPL:ILLUMINA" | samtools view -F 0x04 -q ${QUAL} -bh -U ${OUTPUT_NOMAPPING}${ID}below_q${QUAL}.bam | samtools sort -o ${OUTPUT_MAPPING}${ID}.filt_sort.bam

done
