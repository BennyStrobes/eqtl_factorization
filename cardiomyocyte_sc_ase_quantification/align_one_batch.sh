#!/bin/bash
#SBATCH --time=30:00:00 --mem=20G

batch_name="$1"
bam_input_directory="$2"
vcf_file="$3"
master_barcode_file="$4"
genome_build_file="$5"
gatk_jar="$6"
picard_jar="$7"
ase_read_count_dir="$8"


bam_file=$bam_input_directory"reorderedfinal_"$batch_name".bam"

batch_barcode_file=$ase_read_count_dir$batch_name"_barcodes.txt"
if false; then
python generate_barcodes_for_this_batch.py $batch_name $master_barcode_file $batch_barcode_file
fi


output_prefix=$ase_read_count_dir$batch_name
if false; then
python split_bams_by_barcode.py $bam_file $batch_barcode_file $output_prefix "."
fi


# Loop through cell lines
if false; then
while IFS=$'\t' read -r -a myArray
do
	full_cell_id=${myArray[0]}
	line_id=${myArray[2]}
	if [[ $full_cell_id = *$batch_name* ]]
	then 
		IFS='_' read -r -a array <<< "$full_cell_id"
		cell_id=${array[1]}"_"${array[2]}

		IFS='A' read -r -a line_arr <<< "$line_id"
		short_line_id=${line_arr[1]}
		
		# run ase analysis for this cell line
		bam_file=$ase_read_count_dir$cell_id".bam"
		bam2_file=$ase_read_count_dir$cell_id"_2.bam"
		bam3_file=$ase_read_count_dir$cell_id"_3.bam"
		bam4_file=$ase_read_count_dir$cell_id"_4.bam"
		bam5_file=$ase_read_count_dir$cell_id"_5.bam"
		bam6_file=$ase_read_count_dir$cell_id"_6.bam"
		bam7_file=$ase_read_count_dir$cell_id"_7.bam"
		output_file=$ase_read_count_dir$cell_id"_read_counter_table.txt"

		echo "Starting "$cell_id
		# This is only done so GATK accepts the bams. The group names mean nothing!
		java -jar $picard_jar AddOrReplaceReadGroups \
     		 I=$bam_file  \
     		 O=$bam2_file \
     		 RGID=4 \
     		 RGLB=lib1 \
     		 RGPL=illumina \
     		 RGPU=unit1 \
     		 RGSM=20
	
		samtools sort -o $bam3_file $bam2_file
		samtools index $bam3_file

		samtools calmd -b $bam3_file $genome_build_file > $bam4_file

		samtools sort -o $bam5_file $bam4_file
		samtools index $bam5_file

		# Reordering bam file using picard: reorder sam (ensures correct format for GATK)
		java -jar $picard_jar ReorderSam \
    		I=$bam5_file \
   			O=$bam6_file \
    		R=$genome_build_file \
    		CREATE_INDEX=TRUE 

		metrics_file=$ase_read_count_dir$cell_id"metrics.txt"
		java -jar $picard_jar MarkDuplicates I=$bam6_file METRICS_FILE=$metrics_file O=$bam7_file


		python allelecounter.py --vcf $vcf_file --sample $short_line_id --bam $bam7_file --ref $genome_build_file --min_cov 1 --min_baseq 2 --min_mapq 10 --output $output_file


	fi
done <$master_barcode_file
fi


if false; then
python summarize_batch_ase_calls.py $ase_read_count_dir $batch_name
fi

Rscript visualize_batch_ase_calls.R $ase_read_count_dir $batch_name



























bam_file=$ase_read_count_dir"E1CD2col3_TTTTTTGGAAAT.bam"
bam2_file=$ase_read_count_dir"E1CD2col3_TTTTTTGGAAAT_2.bam"
bam3_file=$ase_read_count_dir"E1CD2col3_TTTTTTGGAAAT_3.bam"
bam4_file=$ase_read_count_dir"E1CD2col3_TTTTTTGGAAAT_4.bam"
bam5_file=$ase_read_count_dir"E1CD2col3_TTTTTTGGAAAT_5.bam"
bam6_file=$ase_read_count_dir"E1CD2col3_TTTTTTGGAAAT_6.bam"
bam7_file=$ase_read_count_dir"E1CD2col3_TTTTTTGGAAAT_7.bam"


output_file=$ase_read_count_dir"E1CD2col3_TTTTTTGGAAAT_gatk_table4.txt"

if false; then
echo "Adding uninformative groupnames to bam file (so GATK accepts the bams)"
#  This is only done so GATK accepts the bams. The group names mean nothing!
java -jar $picard_jar AddOrReplaceReadGroups \
      I=$bam_file  \
      O=$bam2_file \
      RGID=4 \
      RGLB=lib1 \
      RGPL=illumina \
      RGPU=unit1 \
      RGSM=20
fi
if false; then
samtools sort -o $bam3_file $bam2_file
samtools index $bam3_file

samtools calmd -b $bam3_file $genome_build_file > $bam4_file


echo "Reordering bam file using picard: reorder sam (ensures correct format for GATK)"
java -jar $picard_jar ValidateSamFile \
    I=$bam4_file \
    MODE=SUMMARY

samtools sort -o $bam5_file $bam4_file
samtools index $bam5_file

echo "Reordering bam file using picard: reorder sam (ensures correct format for GATK)"
java -jar $picard_jar ReorderSam \
    I=$bam5_file \
    O=$bam6_file \
    R=$genome_build_file \
    CREATE_INDEX=TRUE 
fi

if false; then
python rmdup.py $bam6_file $bam7_file
fi 
metrics_file=$ase_read_count_dir"metrics.txt"
if false; then
java -jar $picard_jar MarkDuplicates I=$bam6_file METRICS_FILE=$metrics_file O=$bam7_file
fi


if false; then
samtools index $bam7_file
fi

echo "Running GATK"
if false; then
java -jar $gatk_jar \
   -R $genome_build_file \
   -T ASEReadCounter \
   -o $output_file \
   -I $bam7_file \
   -sites $vcf_file \
   -U ALLOW_N_CIGAR_READS \
   --outputFormat "TABLE"
fi
metrics_file=$ase_read_count_dir"metrics.txt"
if false; then
java -jar $picard_jar MarkDuplicates I=$bam_file METRICS_FILE=$metrics_file O=$bam2_file
fi

sample_name="18912"
if false; then
python allelecounter.py --vcf $vcf_file --sample $sample_name --bam $bam6_file --ref $genome_build_file --min_cov 1 --min_baseq 2 --min_mapq 10 --output $ase_read_count_dir"temp3"
fi
