#!/bin/bash

SAMTOOLS_DIR="/home-user/software/local/bin"

REF_LEN=3338036
DIR_IN="02_GATK"

BAMS=($(find $DIR_IN -maxdepth 1 -type f -name "*dedup.rnd0.bam"))

for BAM in "${BAMS[@]}"
do
  SUB=$(echo $BAM | sed 's/.*\///g' | sed 's/rnd0.bam$/slm/g')
  OUT=$(echo $BAM | sed 's/bam$/depth/g')

  echo -e "#!/bin/bash -l\n\n#SBATCH -n 1\n" > $SUB
  echo -e "time $SAMTOOLS_DIR/samtools depth $BAM > $OUT\n" >> $SUB
  echo -e "time cnt=\$(awk 'BEGIN{c=0}{if(\$3>=5){c++}}END{print c}' $OUT)" >> $SUB
  echo -e "time let cov=100*\$cnt/$REF_LEN" >> $SUB
  echo -e "echo -e \"good_base\\\tdepth>=5\\\t\$cnt\\\t\$cov\" >> $OUT\n" >> $SUB
  echo -e "time cnt=\$(awk 'BEGIN{c=0}{if(\$3>=10){c++}}END{print c}' $OUT)" >> $SUB
  echo -e "time let cov=100*\$cnt/$REF_LEN" >> $SUB
  echo -e "echo -e \"good_base\\\tdepth>=10\\\t\$cnt\\\t\$cov\" >> $OUT\n" >> $SUB

  chmod +x $SUB
  sbatch $SUB
done

# grep "depth>=10" 02_GATK/*depth | awk '{if($4>=70){print $0}}' |wc -l
# 92
# grep "depth>=10" 02_GATK/*depth | awk '{if($4>=80){print $0}}' |wc -l
# 92
# grep "depth>=10" 02_GATK/*depth | awk '{if($4>=90){print $0}}' |wc -l
# 83

