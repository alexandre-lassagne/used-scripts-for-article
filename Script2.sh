
#!/bin/bash

# Function to display usage message
usage() {
    echo "Usage: $0 
        This script allows to perfom flagstat and quality control mapping
        It will create 2 files in the output_path: 
         - flagstat =>  flagstat.txt
                        depth_a.txt: all positions
                        depth.txt: only positions covered
         - QCResult =>  nbreads_mapped.txt: number of mapped reads on the reference
                        nbreads_tot.txt: number of total reads
                        depth_cover.txt: mean reads detph among positions covered / mean reads depth among all positions / mean reads breath
                        multiqc.html
         
        Default threads = 6
        
        Flags option:
    [-t threads] 
    -i input_mapping 
    -f fastq_path 
    -o output_path
    "
    exit 1
}


# Default value for THREADS
THREADS=6


# Parse command-line arguments
while getopts t:i:f:o:h flag
do
    case "${flag}" in
        t) THREADS=${OPTARG};;
        i) INPUT=${OPTARG};;
        f) PATH_FASTQ=${OPTARG};;
        o) OUTPUT=${OPTARG};;
        h) usage;;
        *) usage;;
    esac
done

# Check if all required parameters are provided
if [ -z "$PATH_FASTQ" ] || [ -z "$OUTPUT" ] || [ -z "$INPUT" ]; then
    echo "Error: Missing required arguments."
    usage
fi

mkdir -p ${OUTPUT}/flagstat
mkdir -p ${OUTPUT}/QCResult


#### Run samtools flagstat and samtools depth on filtered reads
for bam in ${INPUT}/*.filt_sort.bam; do
ID=`basename -s .filt_sort.bam $bam`
samtools flagstat -@ ${THREADS} ${INPUT}/${ID}.filt_sort.bam > ${OUTPUT}/flagstat/${ID}_flagstat.txt
samtools depth -a ${INPUT}/${ID}.filt_sort.bam > ${OUTPUT}/flagstat/${ID}_depth_a.txt
samtools depth ${INPUT}/${ID}.filt_sort.bam > ${OUTPUT}/flagstat/${ID}_depth.txt
done



multiqc ${OUTPUT}/flagstat/*_flagstat.txt  -o ${OUTPUT}/QCResult/

#### Estimate mean read depth
# size_ref= total number of position of ref genome size 
# sum= sum of depth estimates at each position 
# the sum is done only with depth estimated at mapped positions

#### Estimate breath depth = proportion of reads that mapped to the reference genome (_depth_a.txt file shows all positions)
# size_ref= total number of positions of the ref ~ approximation of the the reference genome size
# mapped= number of positions of mapped genome (increment 1 each time there is at least 1 read that mapped at a position)

echo -e "ID\tmean_read_DP\tbreath_DP" > ${OUTPUT}/QCResult/depth_cover.txt
for depth_file in ${OUTPUT}/flagstat/*_depth_a.txt; do
ID=`basename -s _depth_a.txt $depth_file`
mean_read_DPcover=`awk -v OFS='\t' -v spname="$ID" '{size_ref++; if($3>0) sum+=$3}END{print spname, sum/size_ref}' ${OUTPUT}/flagstat/${ID}_depth.txt`
mean_read_DPall=`awk -v OFS='\t' -v spname="$ID" '{size_ref++; if($3>0) sum+=$3}END{print spname, sum/size_ref}' ${OUTPUT}/flagstat/${ID}_depth_a.txt`
breath_DP=`awk -v OFS='\t' -v spname="$ID" '{size_ref++; if($3>0) mapped+=1}END{print(mapped/size_ref)*100}' ${OUTPUT}/flagstat/${ID}_depth_a.txt`
paste <(printf %s "$mean_read_DPcover") <(printf %s "$mean_read_DPall") <(printf %s "$breath_DP")>> ${OUTPUT}/QCResult/depth_cover.txt
done

#### Estimate the proportion of mapped reads to the reference

# 1st way = from samtools flagstat file if you keep unfiltered/on unfiltered reads
# 2nd way = use the number of mapped reads from flagstat filtered file and count the total number of reads in the fastq.gz file
echo -e "ID\tMAPPED_reads"> ${OUTPUT}/QCResult/nbreads_mapped.txt
for flagstat_file in ${OUTPUT}/flagstat/*_flagstat.txt; do
ID=`basename -s _flagstat.txt $flagstat_file`
awk -v spname="$ID" -F "[+]" 'NR == 7 {print spname, $1}' ${OUTPUT}/flagstat/${ID}_flagstat.txt >> ${OUTPUT}/QCResult/nbreads_mapped.txt
done

#in a fastq file, 1read has 4 lines: (nbR1+nbR2)/4 ==> nbR1/2 because nbR1=nbR2
echo -e "ID\tTOTAL_reads"> ${OUTPUT}/QCResult/nbreads_tot.txt
for fastq in ${PATH_FASTQ}*_R1.fastq.gz; do
ID=${fastq/_R1.fastq.gz/}
nbR=`zcat ${ID}_R1.fastq.gz | wc -l`
TOT=$(($nbR/2))
paste <(printf %s "$ID") <(printf %s "$TOT")>>${OUTPUT}/QCResult/nbreads_tot.txt
done

paste ${OUTPUT}/QCResult/nbreads_mapped.txt ${OUTPUT}/QCResult/nbreads_tot.txt > ${OUTPUT}/QCResult/prop_mapped.txt
echo -e "ID\tMAPPED_reads\tID\tTOT\tPROPORTION" > ${OUTPUT}/QCResult/propmapped_final.txt
awk 'NR>1{PROPORTION = 100*$2/$4 ; print $1"\t"$2"\t"$3"\t"$4"\t"PROPORTION}' ${OUTPUT}/QCResult/prop_mapped.txt >> ${OUTPUT}/QCResult/propmapped_final.txt
