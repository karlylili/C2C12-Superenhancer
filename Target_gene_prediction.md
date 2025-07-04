##Target gene prediction
python3 ABC-Enhancer-Gene-Prediction-master/src/makeCandidateRegions.py --bam Starr-seq.bam --chrom_sizes mm10.fai --outDir ./ --nStrongestPeaks 150000 --narrowPeak SE_uniq.narrowPeak 
python3 ABC-Enhancer-Gene-Prediction-master/src/run.neighborhoods.py --candidate_enhancer_regions SE_uniq.narrowPeak.candidateRegions.bed --genes /mm10.gene.bed --ATAC ATAC.bam --expression_table tpm.txt --chrom_sizes mm10.fai --outdir ./ 
python3 ABC-Enhancer-Gene-Prediction-master/src/predict.py --enhancers EnhancerList.txt  --genes GeneList.txt --chrom_sizes mm10.fai --outdir ./ --make_all_putative --threshold 0.02  --HiCdir ./hic/ --hic_resolution 5000
