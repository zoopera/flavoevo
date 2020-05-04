#!/bin/bash -l
#SBATCH -n 1

JAVA_DIR="/home-user/software/jre/bin"
GATK_DIR="/home-user/software/gatk/latest"

DIR_IN="02_GATK"
DIR_OUT="02_GATK"

mkdir -p $DIR_OUT # create the output directory

BOOT=3

# Based on the Long's paper and GATK best practice recommendations
# Float numbers are required
SNP_QUAL=100.0
SNP_MQ=59.0
# SNP_SOR=3.0 # not used by Long
INDEL_QUAL=100.0
INDEL_MQ=59.0
# INDEL_SOR=10.0 # not used by Long

VCFS=($(find $DIR_IN -maxdepth 1 -type f -name "*rnd$BOOT.vcf"))

for VCF in "${VCFS[@]}"
do
  HEAD=$(echo $VCF | sed 's/.*\///g' | sed 's/.vcf$//g')
  SUB=$HEAD".slm"
  echo -e "#!/bin/bash\n" > $SUB
  echo -e "export PATH=$JAVA_DIR:\$PATH\n" >> $SUB

  # filter SNPs
  SNP_IN=$DIR_OUT"/"$HEAD".snp.vcf"
  SNP_OUT=$DIR_OUT"/"$HEAD".snp.flt.vcf"
  echo -e "time $GATK_DIR/gatk VariantFiltration --output $SNP_OUT --variant $SNP_IN --filter-expression \"QUAL <= $SNP_QUAL || MQ <= $SNP_MQ\" --filter-name \"snp_filter\"\n" >> $SUB
  # filter INDELs
  INDEL_IN=$DIR_OUT"/"$HEAD".indel.vcf"
  INDEL_OUT=$DIR_OUT"/"$HEAD".indel.flt.vcf"
  echo -e "time $GATK_DIR/gatk VariantFiltration --output $INDEL_OUT --variant $INDEL_IN --filter-expression \"QUAL <= $INDEL_QUAL || MQ <= $INDEL_MQ\" --filter-name \"indel_filter\"\n" >> $SUB

  chmod +x $SUB
  sbatch $SUB
done

