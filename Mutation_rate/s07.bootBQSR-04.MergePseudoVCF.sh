#!/bin/bash

JAVA_DIR="/home-user/software/jre/bin/"
GATK_DIR="/home-user/software/gatk/latest"

DIR_IN="02_GATK"
DIR_OUT="02_GATK"

mkdir -p $DIR_OUT # create the output directory

# CombineVCFs
BOOT=1
VCFS=($(find $DIR_IN -maxdepth 1 -type f -name "*dedup.rnd$BOOT\.pseudo.vcf"))
LIST=$DIR_OUT"/vcfs.list"

rm -f $LIST
for VCF in "${VCFS[@]}"
do
  echo $VCF >> $LIST
done

HEAD="Combo"
SUB=$HEAD".slm"
CMB=$DIR_OUT"/"$HEAD".rnd"$BOOT".pseudo.vcf"

echo -e "#!/bin/bash -l\n\n#SBATCH -n 1\n" > $SUB
echo -e "export PATH=$JAVA_DIR:\$PATH\n" >> $SUB
echo -e "time $GATK_DIR/gatk MergeVcfs --INPUT $LIST --OUTPUT $CMB" >> $SUB

chmod +x $SUB
./$SUB
