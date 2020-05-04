#!/bin/bash

JAVA_DIR="/home-user/software/jre/bin/"
GATK_DIR="/home-user/software/gatk/latest"

DIR_IN="02_GATK"
DIR_OUT="02_GATK"

mkdir -p $DIR_OUT # create the output directory

# AnalyzeCovariates
BQSRS=($(find $DIR_IN -maxdepth 1 -type f -name "*dedup.rnd0.recal"))

for BQSR in "${BQSRS[@]}"
do
  HEAD=$(echo $BQSR | sed 's/.*\///g' | sed "s/.rnd0.recal$//g")
  SUB=$HEAD".slm"
  BEFORE=$DIR_OUT"/"$HEAD".rnd1.recal"
  AFTER=$DIR_OUT"/"$HEAD".rnd2.recal"
  PLOT=$DIR_OUT"/"$HEAD".recal012.pdf"

  echo -e "#!/bin/bash -l\n\n#SBATCH -n 1\n" > $SUB
  echo -e "export PATH=$JAVA_DIR:\$PATH\n" >> $SUB
  echo -e "time $GATK_DIR/gatk AnalyzeCovariates -bqsr $BQSR --before-report-file $BEFORE --after-report-file $AFTER --plots-report-file $PLOT" >> $SUB

  chmod +x $SUB
  ./$SUB
done

