#!/bin/bash

JAVA_DIR="/home-user/software/jre/bin/"
GATK_DIR="/home-user/software/gatk/latest"

DIR_IN="02_GATK"
DIR_OUT="02_GATK"
REF_GNM="00_ASSEMBLY/Sulfitobacter_sp_EE-36.prokka/Sulfitobacter_sp_EE-36.fna"

BOOT=1
SNP_VCF=$DIR_IN"/Combo.rnd"$BOOT".snp.pseudo.flt.pass.vcf" # produced by s06.bootBQSR-06.VariantFiltration.part2.pl
INDEL_VCF=$DIR_IN"/Combo.rnd"$BOOT".indel.pseudo.flt.pass.vcf"

mkdir -p $DIR_OUT # create the output directory

# BaseRecalibrator
# https://gatkforums.broadinstitute.org/firecloud/discussion/comment/28256#Comment_28256
# To be clear, each iteration of BQSR should be done on the original bam file, not on the output of the previous recalibration.
# The only thing that changes is the set of variants used as known set.
BAMS=($(find $DIR_IN -maxdepth 1 -type f -name "*dedup.rnd$BOOT\.bam"))

for BAM in "${BAMS[@]}"
do
  HEAD=$(echo $BAM | sed 's/.*\///g' | sed "s/.rnd$BOOT\.bam$//g")
  SUB=$HEAD".slm"
  OUT=$DIR_OUT"/"$HEAD".rnd"$BOOT".recal" # recalibrated report

  echo -e "#!/bin/bash -l\n\n#SBATCH -n 1\n" > $SUB
  echo -e "export PATH=$JAVA_DIR:\$PATH\n" >> $SUB
  echo -e "time $GATK_DIR/gatk BaseRecalibrator --input $BAM --output $OUT --reference $REF_GNM --known-sites $SNP_VCF --known-sites $INDEL_VCF" >> $SUB

  chmod +x $SUB
  sbatch $SUB
done

