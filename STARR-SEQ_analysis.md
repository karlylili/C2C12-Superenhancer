##STARR-SEQ alignment
bowtie2_index=mm10
less sample.lst |while read line 
do 
bowtie2  -p 30  --very-sensitive -X 2000 -x  $bowtie2_index -1 ../01_clean/clean_data/${line}_1_val_1.fq.gz -2 ../01_clean/clean_data/${line}_2_val_2.fq.gz \|samtools sort  -O bam  -@ 5 -o ${line}.raw.bam 
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
nohup macs2 callpeak -t ${line} -c ORI-mouse-input3.bam  -g 2725521370 -p 0.01 --nomodel --shift 0 --extsize 150 -B --SPMR --keep-dup 1 -n ${line} --outdir ./peak2 &
done 


#summary reads count
multiBamSummary BED-file --BED WT.bed --bamfiles WT_all_merge.bam ../ORI-mouse-input3.bam -out WT_inputCounts.npz --outRawCounts WT_inputCounts &

##activity calculation using R
data2<-read.table("WT_inputCounts",header=T,sep="\t")
data2$fc_a <- data2$output_a.bam/"fragment numbers of a"*1000000/data$input.bam/1000000

##activity hist plot
ggplot(data2, aes(x = mean_FC)) +geom_histogram(binwidth = 0.1, fill = "lightblue", colour = "black")+
geom_vline(aes(xintercept=1), colour="#BB0000", linetype="dashed")+xlim(0,6)+
xlab("Silencer activity")+ theme(axis.text.x=element_text(angle=15,hjust = 1,colour="black",family="Times",size=20),
        axis.text.y=element_text(family="Times",size=16,face="plain"), 
        axis.title.y=element_text(family="Times",size = 20,face="plain"), 
        axis.title.x=element_text(family="Times",size = 20,face="plain"), 
       legend.text=element_text(face="italic", family="Times", colour="black",  
                                 size=16),
        legend.title=element_text(face="italic", family="Times", colour="black", 
                                  size=18))+theme_bw()
##motif enrichment
findMotifsGenome.pl peak.file genome.fa output_motif -len 8,10,12 -size -100,100

##signal distribution
computeMatrix reference-point -S sample.bw -R day0_A1 day0_A2 --skipZeros -omatrix.gz --referencePoint center -a 1000 -b 1000


##enhancer annotation using ChIPseeker
library(ChIPseeker)
library(clusterProfiler)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(GenomeInfoDb)
peak <- readPeakFile("enhancer.bed")
peakAnno <- annotatePeak( peak , tssRegion = c(-3000, 3000), TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene, annoDb = "org.Mm.eg.db")
write.table(as.data.frame(peakAnno),"peak_peak.annotation.tsv", sep="\t",row.names = F,quote = F)
pdf(file = "peak_plot.pdf")
plotAnnoPie(peakAnno)
dev.off()
