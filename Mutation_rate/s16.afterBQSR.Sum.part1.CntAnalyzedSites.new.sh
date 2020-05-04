#!/bin/bash

dir="02_GATK"
cnts=($(find $dir -maxdepth 1 -type f -name "*.clean.sort.dedup.rnd3.cnt"))

for cnt in "${cnts[@]}"
do
  sub=$(echo $cnt | sed 's/.*\///g' | sed 's/cnt$/slm/g')
  echo -e "#!/bin/bash\n" > $sub
  echo -e "time ./s16.afterBQSR.Sum.part1.CntAnalyzedSites.new.pl $cnt" >> $sub
  chmod +x $sub
  sbatch $sub
done
