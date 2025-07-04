## Super enhancer identification
python ABC-Enhancer-Gene-Prediction-master/src/makeCandidateRegions.py --narrowPeak peaks.narrowPeak.sorted --bam Starr-seq.bam --chrom_sizes mm10.fai --outDir ./ --nStrongestPeaks 150000 
python ABC-Enhancer-Gene-Prediction-master/src/run.neighborhoods.py --candidate_enhancer_regions  peaks.narrowPeak.sorted.candidateRegions.bed --genes mm10.gene.bed --ATAC ATAC.bam  --expression_table tpm.txt --chrom_sizes mm10.fai --outDir ./ 
python ABC-Enhancer-Gene-Prediction-master/src/predict.py --enhancers EnhancerList.txt --genes GeneList.txt --chrom_sizes mm10.fai --outdir ./ --make_all_putative --threshold 0.02 --score_column powerlaw.Score 
python ABC-Enhancer-Gene-Prediction-master/src/predict.py --enhancers EnhancerList.txt --genes GeneList.txt --chrom_sizes mm10.fai --outdir ./ --make_all_putative --threshold 0.02  --HiCdir ./hic/ --hic_resolution 5000

