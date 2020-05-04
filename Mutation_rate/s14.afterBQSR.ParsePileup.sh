#!/bin/bash

dir="02_GATK"
pups=($(find $dir -maxdepth 1 -type f -name "*pileup"))

for pup in "${pups[@]}"
do
  sub=$(echo $pup | sed 's/.*\///g' | sed 's/pileup$/slm/g')

  echo -e "#!/bin/bash\n" > $sub
  echo -e "time ./s14.afterBQSR.ParsePileup.pl $pup\n" >> $sub

  chmod +x $sub
  sbatch $sub
done
