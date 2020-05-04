#!/bin/bash

REF="00_ASSEMBLY/Sulfitobacter_sp_EE-36.prokka/Sulfitobacter_sp_EE-36.fna"
DIR_IN="00_READS"
DIR_OUT="01_BWA"
FRS=($(find $DIR_IN -maxdepth 1 -type f -name "*_R1_001.trm.pe.fq.gz"))

mkdir -p $DIR_OUT

# Index database sequences in the FASTA format.
time /home-user/software/bwa/latest/bwa index $REF

# BWA mem
for FR in "${FRS[@]}"
do
  RR=$(echo $FR | sed 's/_R1_001.trm.pe.fq.gz/_R2_001.trm.pe.fq.gz/g')
  SM=$(echo $FR | sed 's/.*\///g' | sed 's/_R1_001.trm.pe.fq.gz/.sam/g')
  SUB=$(echo $FR | sed 's/.*\///g' | sed 's/_R1_001.trm.pe.fq.gz/.slm/g')

  echo -e "#!/bin/bash -l\n\n#SBATCH -n 1\n" > $SUB
  echo -e "time /home-user/software/bwa/latest/bwa mem $REF $FR $RR > $DIR_OUT/$SM" >> $SUB
  chmod +x $SUB
  sbatch $SUB
done
