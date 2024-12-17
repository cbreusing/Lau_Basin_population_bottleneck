!/bin/bash
SBATCH -J Usearch
SBATCH -t 100:00:00
SBATCH -n 32
SBATCH -N 1
SBATCH -p uri-cpu
SBATCH --mem=150g
SBATCH -o Usearch.out
SBATCH -e Usearch.err

module load seqtk/1.3
module load blast-plus/2.13.0+py3.8.12

mkdir usearch_clustering
cd usearch_clustering

eval "$(conda shell.bash hook)"
conda activate vsearch
usearch11 -fastq_mergepairs ../*_paired_R1.fq -fastqout Lau_16S_merged.fq -fastq_minmergelen 230 -fastq_maxmergelen 270 -fastq_maxdiffs 10 -fastq_pctid 80 -relabel @
vsearch --orient Lau_16S_merged.fq --db /project/pi_rbeinart_uri_edu/Databases/SILVA_138.1_SSURef_NR99_tax_silva_trunc.fasta --fastqout Lau_16S_oriented.fq --notmatched Lau_16S_notmatched.fq --notrunclabels
cat Lau_16S_notmatched.fq >> Lau_16S_oriented.fq
vsearch --fastq_filter Lau_16S_oriented.fq --fastq_maxee 1 --fastaout Lau_16S_filtered.fa --fastq_qmax 42 --fastq_minlen 230 --fastq_truncqual 20
vsearch --derep_fulllength Lau_16S_filtered.fa --output Lau_16S_uniques.fa --sizeout --relabel Uniq
usearch11 -unoise3 Lau_16S_uniques.fa -zotus Lau_16S_zotus.fa -tabbedout unoise3.txt
seqtk seq -a Lau_16S_oriented.fq | sed -re "s/>([A-Za-z0-9]+.*)(\.[0-9]+)/>\1\2;sample=\1/g" > Lau_16S_oriented.fa
vsearch --search_exact Lau_16S_oriented.fa --db Lau_16S_zotus.fa --otutabout Lau_16S_exact_zotus.txt --strand plus --threads ${SLURM_CPUS_ON_NODE}
conda deactivate

conda activate qiime2-amplicon-2024.2

qiime tools import --input-path Lau_16S_zotus.fa --output-path Lau_16S_seqs.qza --type 'FeatureData[Sequence]'

qiime feature-classifier extract-reads --i-sequences /project/pi_rbeinart_uri_edu/Databases/silva-138-99-ref-seqs.qza --p-f-primer GTGYCAGCMGCCGCGGTAA --p-r-primer GGACTACNVGGGTWTCTAAT --p-min-length 230 --p-max-length 270 --o-reads silva-138-99-ref-seqs.qza
qiime feature-classifier fit-classifier-naive-bayes --i-reference-reads silva-138-99-ref-seqs.qza --i-reference-taxonomy /project/pi_rbeinart_uri_edu/Databases/silva-138-99-tax.qza --o-classifier silva-138-99-515F-806R-uniform-classifier.qza
qiime clawback assemble-weights-from-Qiita --i-classifier silva-138-99-515F-806R-uniform-classifier.qza --i-reference-taxonomy /project/pi_rbeinart_uri_edu/Databases/silva-138-99-tax.qza --i-reference-sequences silva-138-99-ref-seqs.qza --p-metadata-key empo_3 --p-metadata-value "Surface (saline)" --p-metadata-value "Water (saline)" --p-context Deblur_2021.09-Illumina-16S-V4-150nt-ac8c0b --o-class-weight marine-habitat-weights.qza
qiime feature-classifier fit-classifier-naive-bayes --i-reference-reads silva-138-99-ref-seqs.qza --i-reference-taxonomy /project/pi_rbeinart_uri_edu/Databases/silva-138-99-tax.qza --i-class-weight marine-habitat-weights.qza --o-classifier silva-138-99-515F-806R-weighted-classifier.qza 

qiime feature-classifier classify-sklearn --i-classifier silva-138-99-515F-806R-weighted-classifier.qza --i-reads Lau_16S_seqs.qza --o-classification Lau_16S_taxonomy.qza --p-n-jobs ${SLURM_CPUS_ON_NODE}

qiime tools export --input-path Lau_16S_taxonomy.qza --output-path qiime2_results

conda deactivate

mv qiime2_results/taxonomy.tsv Lau_16S_taxonomy.tsv
rm -r qiime2_results

grep "Gammaproteobacteria" Lau_16S_taxonomy.tsv > Lau_16S_taxonomy_filtered.tsv
grep "Campylobacteria" Lau_16S_taxonomy.tsv >> Lau_16S_taxonomy_filtered.tsv 

perl -anle 'print $F[0]' Lau_16S_taxonomy_filtered.tsv > potential_symbiont_ASVs.txt

seqtk subseq Lau_16S_zotus.fa potential_symbiont_ASVs.txt > Lau_16S_zotus_filtered.fa

blastn -query Lau_16S_zotus_filtered.fa -db Symbiont_refseqs.fasta -out blast.txt -evalue 1e-20 -num_threads ${SLURM_CPUS_ON_NODE} -max_target_seqs 1 -outfmt "6 std stitle"
cat blast.txt | sort -k1,1 -k12,12nr -k11,11n | sort -u -k1,1 > blast.topHit.txt

perl -anle 'print $_ if $F[2]>99' blast.topHit.txt > blast.symbionts.txt

perl -anle 'print $F[0]' blast.symbionts.txt > symbiont_ASVs.txt

seqtk subseq Lau_16S_zotus_filtered.fa symbiont_ASVs.txt > Lau_16S_symbionts.fa
grep -w -f symbiont_ASVs.txt Lau_16S_exact_zotus.txt > Lau_16S_exact_symbionts.txt
grep -w -f symbiont_ASVs.txt Lau_16S_taxonomy_filtered.tsv > Lau_16S_taxonomy_symbionts.tsv
