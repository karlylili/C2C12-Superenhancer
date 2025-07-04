##HiC alignment 
bin/HiC-Pro -i hic_data/ -o result -c config-hicpro.txt

##file format transfer
hic2cool convert -r 25000 day0.allValidPairs.hic day0_25k.cool 
hicConvertFormat  -m day0_50k.cool --inputFormat cool --outputFormat h5 --resolutions 50000 --outFileName day0_50k.h5 

##find TAD 
hicFindTADs -m day4_50k.h5 --outPrefix day4_50k --correctForMultipleTesting fdr -p 10

##differenial TAD
hicDifferentialTAD -tm ../day2/day2_50k.h5 -cm ../day0/day0_50k.h5 -td ./day2/day2_50k_domains.bed -o diff_tad_day0-day2 -t 10 

##find loop
java -Xmx30g -jar juicer_tools.jar hiccups --cpu --threads 8 -r 30000 -k KR day0.allValidPairs.hic ./day0

##differential Loop
java -jar /data2/liguoli/software/juicer-1.6/scripts/common/juicer_tools.jar hiccupsdiff ../../day0/day0.allValidPairs.hic day0_merged_loo.bedpe day2_merged_loops.bedpe ./day0_day2_diff/ --cpu

##AB compartment
makeTagDirectory test -format HiCsummary test.homer
runHiCpca.pl test-500k test -cpu 32 -res 500000 -genome mm10.fasta -pc 1
