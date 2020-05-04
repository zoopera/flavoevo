#!/bin/bash -l
#SBATCH -n 1

JAVA_DIR="/home-user/software/jre/bin/"
GATK_DIR="/home-user/software/gatk/latest"

DIR_IN="02_GATK"
DIR_OUT="02_GATK"

mkdir -p $DIR_OUT # create the output directory

# CombineVCFs
BOOT=2
SNP_IN=$DIR_IN"/Combo.rnd"$BOOT".snp.pseudo.vcf"
SNP_OUT=$DIR_OUT"/Combo.rnd"$BOOT".snp.pseudo.flt.vcf"
INDEL_IN=$DIR_IN"/Combo.rnd"$BOOT".indel.pseudo.vcf"
INDEL_OUT=$DIR_OUT"/Combo.rnd"$BOOT".indel.pseudo.flt.vcf"

# Based on the results of s07.bootBQSR-05.PseudoVCFStats.pl
# Float numbers are required
SNP_QUAL=100.0
SNP_MQ=40.0
SNP_SOR=3.0
INDEL_QUAL=100.0
INDEL_MQ=40.0
INDEL_SOR=10.0

time $GATK_DIR/gatk VariantFiltration --output $SNP_OUT --variant $SNP_IN --filter-expression "QUAL < $SNP_QUAL || MQ < $SNP_MQ || SOR > $SNP_SOR" --filter-name "snp_filter"
time $GATK_DIR/gatk VariantFiltration --output $INDEL_OUT --variant $INDEL_IN --filter-expression "QUAL < $INDEL_QUAL || MQ < $INDEL_MQ || SOR > $INDEL_SOR" --filter-name "indel_filter"

