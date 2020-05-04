#!/bin/bash

DIR_IN="/home-db/luolab/mutation_accumulation_db/2018_sulfitobacter_sp_ee-36/20180723_unh_sulfitobacter_sp_ee-36.re-sequencing/" # change the dir
DIR_OUT="00_READS"

mkdir -p $DIR_OUT

FQS=($(find $DIR_IN -type f -name "*.fastq.gz"))

for FQ in "${FQS[@]}"
do
	SUB=$(echo $FQ|sed 's/.*\///g'|sed 's/gz$/slm/g')
	echo -e "#!/bin/bash -l\n#SBATCH -n 1\n" > $SUB
	echo -e "export PATH=/home-user/software/jre/bin/:\$PATH\n" >> $SUB
	echo -e "time /home-user/software/fastqc/latest/fastqc -o $DIR_OUT $FQ" >> $SUB
	chmod +x $SUB
	sbatch --exclude=cl004 $SUB # node cl004 has problem with java
done
