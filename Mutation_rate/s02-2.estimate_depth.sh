#!/bin/bash

DIR_IN="00_READS"

GSIZE=3338036

FREADS=($(find $DIR_IN -maxdepth 1 -type f -name "*R1_001.trm.pe.fq.gz"))

for FREAD in "${FREADS[@]}"
do
	TOTAL=0
	CNT=$(zcat $FREAD|paste - - - -|cut -f2|awk 'BEGIN{c=0}{c+=length($0)}END{print c}')
	let TOTAL=$TOTAL+$CNT
	r=$(echo $FREAD|sed 's/R1_001.trm.pe.fq.gz/R2_001.trm.pe.fq.gz/g')
	CNT=$(zcat $r|paste - - - -|cut -f2|awk 'BEGIN{c=0}{c+=length($0)}END{print c}')
	let TOTAL=$TOTAL+$CNT
	let TOTAL=$TOTAL/$GSIZE
	echo -e "$FREAD\t$TOTAL"
done
