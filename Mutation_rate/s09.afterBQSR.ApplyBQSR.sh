#!/bin/bash

JAVA_DIR="/home-user/software/jre/bin/"
GATK_DIR="/home-user/software/gatk/latest"
SAMTOOLS_DIR="/home-user/software/local/bin/"

DIR_IN="02_GATK"
DIR_OUT="02_GATK"

BOOT=2
let NEXT=$BOOT+1

mkdir -p $DIR_OUT # create the output directory

# BaseRecalibrator
BAMS=($(find $DIR_IN -maxdepth 1 -type f -name "*dedup.rnd$BOOT\.bam"))

for BAM in "${BAMS[@]}"
do
  HEAD=$(echo $BAM | sed 's/.*\///g' | sed "s/.rnd$BOOT\.bam$//g")
  SUB=$HEAD".slm"
  RECALB_TAB=$DIR_OUT"/"$HEAD".rnd"$BOOT".recal"
  RECALB_BAM=$DIR_OUT"/"$HEAD".rnd"$NEXT".bam"

  echo -e "#!/bin/bash -l\n\n#SBATCH -n 1\n" > $SUB
  echo -e "export PATH=$JAVA_DIR:\$PATH\n" >> $SUB
  echo -e "time $GATK_DIR/gatk ApplyBQSR --input $BAM --output $RECALB_BAM --bqsr-recal-file $RECALB_TAB" >> $SUB

  chmod +x $SUB
  sbatch $SUB
done

