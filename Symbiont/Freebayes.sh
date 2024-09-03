#!/bin/bash
#SBATCH -J Freebayes
#SBATCH -t 100:00:00
#SBATCH -n 24
#SBATCH -N 1
#SBATCH -p batch
#SBATCH --mem=100g
#SBATCH -o Freebayes.out
#SBATCH -e Freebayes.err
#SBATCH --account=epscor-condo

module load samtools/1.16.1-txuglks 
module load bowtie2/2.4.2-xdquyzq
module load gatk/4.3.0.0-234wqft 
module load picard/2.26.2-qabtyqy
module load r/4.2.2-z6qdiis
module load vcftools/0.1.14-syssqsi
module load bcftools/1.16-ewu6fpe
module load htslib/1.17-zxcat2k

bowtie2-build Epsilon_pangenome.fasta Epsilon_pangenome.fasta
samtools faidx Epsilon_pangenome.fasta
 
for file in `cat filelist.txt`
do
Map reads to reference
bowtie2 --very-sensitive-local -p ${SLURM_CPUS_ON_NODE} -x Epsilon_pangenome.fasta -1 ../../../${file}_R1_clean.fastq -2 ../../../${file}_R2_clean.fastq | samtools view -bS -h -@ ${SLURM_CPUS_ON_NODE} - | samtools sort -@ ${SLURM_CPUS_ON_NODE} - > ${file}.Epsilon.sorted.bam
Mark and remove duplicates
picard MarkDuplicates I=${file}.Epsilon.sorted.bam O=${file}.Epsilon.dedup.bam M=${file}.Epsilon.metrics.txt REMOVE_DUPLICATES=true
samtools index ${file}.Epsilon.dedup.bam
Indel realignment
eval "$(conda shell.bash hook)"
conda activate lofreq
lofreq viterbi -f Epsilon_pangenome.fasta -k ${file}.Epsilon.dedup.bam | samtools sort -@ ${SLURM_CPUS_ON_NODE} - > ${file}.Epsilon.realigned.bam
Base recalibration
lofreq indelqual -f Epsilon_pangenome.fasta --dindel -o ${file}.Epsilon.indelqual.bam ${file}.Epsilon.realigned.bam
samtools index ${file}.Epsilon.indelqual.bam
samtools view -bS -h -F4 ${file}.Epsilon.indelqual.bam > ${file}.Epsilon.filt.dedup.bam
Check number of alignments for each BAM file, select the lowest number of alignments for downsampling
samtools view -c ${file}.Epsilon.filt.dedup.bam >> alignment_counts.txt
done

rm *.Epsilon.sorted.bam
rm *.Epsilon.dedup.bam*
rm *.Epsilon.realigned.bam
rm *.Epsilon.indelqual.bam*

num=947663

for file in `cat filelist.txt`
do
count=`samtools view -c ${file}.Epsilon.filt.dedup.bam`
if [ $num -le $count ]
   then
    frac=`bc -l <<< $num/$count`
    samtools view -@ ${SLURM_CPUS_ON_NODE} -h -bs $frac ${file}.Epsilon.filt.dedup.bam > ${file}.Epsilon.subsampled.bam
fi
samtools index ${file}.Epsilon.subsampled.bam
picard AddOrReplaceReadGroups I=${file}.Epsilon.subsampled.bam O=${file}.Epsilon.subsampled.RG.bam RGID=${file} RGLB=LIB_${file} RGPL=ILLUMINA RGPU=FLOWCELL1 RGSM=${file} VALIDATION_STRINGENCY=SILENT
samtools index ${file}.Epsilon.subsampled.RG.bam
bam trimBam ${file}.Epsilon.subsampled.RG.bam - -c | samtools sort -@ ${SLURM_CPUS_ON_NODE} -n - > ${file}.Epsilon.softclipped.bam
samtools fixmate ${file}.Epsilon.softclipped.bam - | samtools sort -@ ${SLURM_CPUS_ON_NODE} - > ${file}.Epsilon.softclipped.sorted.bam
done

rm *.Epsilon.subsampled.bam*
rm *.Epsilon.subsampled.RG.bam*
rm *.Epsilon.softclipped.bam

ls *.Epsilon.softclipped.sorted.bam > bam.list

freebayes -f Epsilon_pangenome.fasta -L bam.list -v Epsilon.Freebayes.vcf -F 0.01 -C 1 -p 1 --pooled-continuous -g 1000 -m 30 -q 20 --min-coverage 10 --haplotype-length 0 --report-monomorphic

cat Epsilon.Freebayes.vcf | vcf-sort -c > Epsilon.Freebayes.sorted.vcf
bgzip -c Epsilon.Freebayes.sorted.vcf > Epsilon.Freebayes.vcf.gz
tabix -p vcf Epsilon.Freebayes.vcf.gz
bcftools filter -g 5 -i 'REF!="N" && SRP > 5 && SAP > 5 && EPP > 5 && QUAL > 20 && INFO/DP > 10' -o Epsilon.Freebayes.filtered.vcf -O v Epsilon.Freebayes.vcf.gz
vcftools --vcf Epsilon.Freebayes.filtered.vcf --max-missing 0.75 --recode --recode-INFO-all --out Epsilon.Freebayes
vcftools --vcf Epsilon.Freebayes.recode.vcf --missing-indv --out Epsilon.Freebayes.recode
vcftools --vcf Epsilon.Freebayes.recode.vcf --site-mean-depth --out Epsilon.Freebayes.recode
vcftools --vcf Epsilon.Freebayes.recode.vcf --max-meanDP 20 --recode --recode-INFO-all --out Epsilon.Freebayes.FINAL
vcftools --vcf Epsilon.Freebayes.FINAL.recode.vcf --extract-FORMAT-info GT --out Epsilon.Freebayes.FINAL
gatk VariantsToTable -V Epsilon.Freebayes.FINAL.recode.vcf -O Epsilon.Freebayes.FINAL.RD.FORMAT -F CHROM -F POS -F REF -ASGF RO
gatk VariantsToTable -V Epsilon.Freebayes.FINAL.recode.vcf -O Epsilon.Freebayes.FINAL.AD.FORMAT -F CHROM -F POS -F ALT -ASGF AO --split-multi-allelic
