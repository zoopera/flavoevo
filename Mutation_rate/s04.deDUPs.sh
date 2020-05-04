#!/bin/bash

JAVA_DIR="/home-user/software/jre/bin/"
GATK_DIR="/home-user/software/gatk/latest"
SAMTOOLS_DIR="/home-user/software/local/bin/"

DIR_IN="01_BWA"
DIR_OUT="02_GATK"

BOOT=0

mkdir -p $DIR_OUT # create the output directory

SAMS=($(find $DIR_IN -maxdepth 1 -type f -name "*sam"))

for SAM in "${SAMS[@]}"
do
  HEAD=$(echo $SAM | sed 's/.*\///g' | sed 's/.sam$//g')
  SUB=$HEAD".slm"
  SAM_RD=$DIR_OUT"/"$HEAD".rhead.sam"
  SAM_CLEAN=$DIR_OUT"/"$HEAD".clean.sam"
  BAM_SORT=$DIR_OUT"/"$HEAD".clean.sort.bam"
  BAM_DEDUP=$DIR_OUT"/"$HEAD".clean.sort.dedup.rnd"$BOOT".bam"
  MTX_DEDUP=$DIR_OUT"/"$HEAD".clean.sort.dedup.rnd"$BOOT".mtx"

  # read group
  RGLB=$(echo $HEAD | sed 's/_.*//g' | sed 's/^/lib/g') # Read-Group library (= sample name)
  RGPL="illumina" # Read-Group platform
  RGPU=$(echo $HEAD | sed 's/_L00.*//g' | sed 's/.*_//g') # Read-Group platform unit
  RGSM=$(echo $HEAD | sed 's/_.*//g' | sed 's/^/smp/g') # Read-Group sample name

  echo -e "#!/bin/bash -l\n\n#SBATCH -n 1\n" > $SUB
  echo -e "export PATH=$JAVA_DIR:\$PATH\n" >> $SUB
  echo -e "time $GATK_DIR/gatk AddOrReplaceReadGroups --INPUT $SAM --OUTPUT $SAM_RD --RGLB=$RGLB --RGPL=$RGPL --RGPU=$RGPU --RGSM=$RGSM" >> $SUB
  echo -e "time $GATK_DIR/gatk CleanSam --INPUT $SAM_RD --OUTPUT $SAM_CLEAN" >> $SUB
  echo -e "time $GATK_DIR/gatk SortSam --INPUT $SAM_CLEAN --OUTPUT $BAM_SORT --SORT_ORDER coordinate" >> $SUB
  echo -e "time $GATK_DIR/gatk MarkDuplicates --INPUT $BAM_SORT --OUTPUT $BAM_DEDUP --METRICS_FILE $MTX_DEDUP --REMOVE_DUPLICATES true\n" >> $SUB
  echo -e "time $SAMTOOLS_DIR/samtools index $BAM_DEDUP" >> $SUB
  chmod +x $SUB
  sbatch $SUB
done

