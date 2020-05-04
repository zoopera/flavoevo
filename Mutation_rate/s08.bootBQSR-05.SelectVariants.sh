#!/bin/bash

JAVA_DIR="/home-user/software/jre/bin/"
GATK_DIR="/home-user/software/gatk/latest"

export PATH=$JAVA_DIR:$GATK_DIR:$PATH

DIR_IN="02_GATK"
DIR_OUT="02_GATK"

mkdir -p $DIR_OUT # create the output directory

BOOT=2
REF_GNM="00_ASSEMBLY/Sulfitobacter_sp_EE-36.prokka/Sulfitobacter_sp_EE-36.fna"
VCF_IN=$DIR_IN"/Combo.rnd"$BOOT".pseudo.vcf"

# select SNPs
SNP_OUT=$DIR_OUT"/Combo.rnd"$BOOT".snp.pseudo.vcf"
gatk SelectVariants --variant $VCF_IN --output $SNP_OUT --select-type-to-include SNP --reference $REF_GNM

# select INDELs
INDEL_OUT=$DIR_OUT"/Combo.rnd"$BOOT".indel.pseudo.vcf"
gatk SelectVariants --variant $VCF_IN --output $INDEL_OUT --select-type-to-include INDEL --reference $REF_GNM


