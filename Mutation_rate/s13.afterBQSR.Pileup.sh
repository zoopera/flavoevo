#!/bin/bash -l
#SBATCH -n 1

JAVA_DIR="/home-user/software/jre/bin"
GATK_DIR="/home-user/software/gatk/latest"
SAMTOOLS_DIR="/home-user/software/local/bin"

DIR_IN="02_GATK"
DIR_OUT="02_GATK"
REF_GNM="00_ASSEMBLY/Sulfitobacter_sp_EE-36.prokka/Sulfitobacter_sp_EE-36.fna" # use the updated reference genome

mkdir -p $DIR_OUT # create the output directory

BOOT=3

BAMS=($(find $DIR_IN -maxdepth 1 -type f -name "*rnd$BOOT.bam"))

for BAM in "${BAMS[@]}"
do
  HEAD=$(echo $BAM | sed 's/.*\///g' | sed 's/.bam$//g')
  SUB=$HEAD".slm"
  PILEUP=$DIR_OUT"/"$HEAD".pileup"

  echo -e "#!/bin/bash\n" > $SUB
  echo -e "export PATH=$JAVA_DIR:\$PATH\n" >> $SUB
  echo -e "time $SAMTOOLS_DIR/samtools mpileup --fasta-ref $REF_GNM -a -d 10000000 --min-MQ 20 --min-BQ 20 --output $PILEUP $BAM\n" >> $SUB

  chmod +x $SUB
  sbatch $SUB
done
