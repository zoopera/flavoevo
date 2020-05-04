#!/bin/bash
#useful when adapter sequences are provided by the platform
#DIR_IN="/home-db/luolab/mutation_accumulation_db/2018_sulfitobacter_sp_ee-36/20180723_unh_sulfitobacter_sp_ee-36.re-sequencing/"
DIR_IN="/home-user/xjwang/Prochlorococcus_MA/00_READS/raw_data"
DIR_OUT="00_READS"

mkdir -p $DIR_OUT

FRS=($(find $DIR_IN -type f -name "M*_1.fq.gz"))

for FR in "${FRS[@]}"
do
	SUB=$(echo $FR | sed 's/.*\///g' | sed 's/_1.*/.slm/g')
	RR=$(echo $FR | sed 's/_1/_2/g')
	ADPT=$(echo $FR | sed 's/.*\///g' | sed 's/_1.fastq.gz/.adapter/g')

	echo -e "#!/bin/bash -l\n#SBATCH -n 1\n" > $SUB
	echo -e "export PATH=/home-user/software/jre/bin/:\$PATH\n" >> $SUB
	echo -e "time /home-user/software/bbmap/bbmerge.sh in1=$FR in2=$RR outa=$DIR_OUT"/"$ADPT reads=-1\n" >> $SUB
	chmod +x $SUB
	sbatch $SUB
done

