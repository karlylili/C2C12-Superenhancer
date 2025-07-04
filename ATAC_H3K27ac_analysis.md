##alignment
bowtie2_index=mm10
less sample.lst |while read line 
do 
bowtie2  -p 30  --very-sensitive -X 2000 -x  $bowtie2_index -1 ../../01_clean/clean_data/${line}_1_val_1.fq.gz -2 ../../01_clean/clean_data/${line}_2_val_2.fq.gz \|samtools sort  -O bam  -@ 5 -o ${line}.raw.bam 
samtools index ${line}.raw.bam 
bedtools bamtobed -i ${line}.raw.bam  \> ${line}.raw.bed
samtools flagstat ${line}.raw.bam  \> ${line}.raw.stat


sambamba markdup --overflow-list-size 600000  --tmpdir='./'  -r ${line}.raw.bam  ${line}.rmdup.bam
samtools index   ${line}.rmdup.bam

samtools flagstat  ${line}.rmdup.bam \> ${line}.rmdup.stat
samtools view  -h  -f 2 -q 30    ${line}.rmdup.bam   \|grep -v chrM \|samtools sort  -O bam  -@ 5 -o  ${line}.last.bam
samtools index   ${line}.last.bam
samtools flagstat  ${line}.last.bam \> ${line}.last.stat
bamCoverage --bam input/${line}.last.bam -o bw/${line}.bw --binSize 10 --extendReads --normalizeUsing RPKM
bamCoverage --bam ${line}.last.bam -o ${line}.bw --binSize 10 --extendReads --normalizeUsing RPKM >> bw.exe
done 


##call peak
macs2 callpeak -t ../01_align/${line} -g  2725521370 -p 0.01 --nomodel --shift -75 --extsize 150 -B --SPMR --keep-dup all --call-summits -n ${line} --outdir ./

##signal distribution
computeMatrix reference-point -S sample.bw -R day0_A1 day0_A2 --skipZeros -omatrix.gz --referencePoint center -a 1000 -b 1000

## motif enrichment 
findMotifsGenome.pl peak.file genome.fa output_motif -len 8,10,12 -size -100,100

##
