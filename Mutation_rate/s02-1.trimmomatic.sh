#!/bin/bash

JAVA_PATH="/home-user/software/jre/bin/"
TRM_PATH="/home-user/software/trimmomatic/v0.36/Trimmomatic-0.36" # where trimmomatic is
DIR_IN="/home-db/luolab/mutation_accumulation_db/2018_sulfitobacter_sp_ee-36/20180723_unh_sulfitobacter_sp_ee-36.re-sequencing/" # change the DIR_IN
DIR_OUT="00_READS"
FREADS=($(find $DIR_IN -type f -name "*R1_001.fastq.gz"))

mkdir -p $DIR_OUT

for FREAD in "${FREADS[@]}"
do
	SUB=$(echo $FREAD | sed 's/.*\///g' | sed 's/_R1_001.fastq.gz$/.slm/g')
	LOG=$(echo $FREAD | sed 's/.*\///g' | sed 's/_R1_001.fastq.gz$/.log/g')
	RREAD=$(echo $FREAD | sed 's/R1_001/R2_001/g')
	ADPT=$(echo $FREAD | sed 's/.*\///g' | sed 's/_R1_001.fastq.gz/.adapter/g')
	FP_OUT=$(echo $FREAD | sed 's/.*\///g' | sed 's/fastq.gz/trm.pe.fq.gz/g') # output_forward_paired.fq.gz
	FS_OUT=$(echo $FREAD | sed 's/.*\///g' | sed 's/fastq.gz/trm.se.fq.gz/g') # output_forward_unpaired.fq.gz 
	RP_OUT=$(echo $RREAD | sed 's/.*\///g' | sed 's/fastq.gz/trm.pe.fq.gz/g') # output_reverse_paired.fq.gz
	RS_OUT=$(echo $RREAD | sed 's/.*\///g' | sed 's/fastq.gz/trm.se.fq.gz/g') # output_reverse_unpaired.fq.gz

	echo -e "#!/bin/bash -l\n#SBATCH -n 1\n" > $SUB
	echo -e "time $JAVA_PATH/java -jar $TRM_PATH/trimmomatic-0.36.jar PE -threads 1 -phred33 -trimlog $LOG $FREAD $RREAD $DIR_OUT/$FP_OUT $DIR_OUT/$FS_OUT $DIR_OUT/$RP_OUT $DIR_OUT/$RS_OUT ILLUMINACLIP:$DIR_OUT/$ADPT:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MAXINFO:40:0.8 MINLEN:40" >> $SUB
	chmod +x $SUB
	sbatch $SUB
done

