#!/bin/bash
#SBATCH -J GenomeWiseCalculations
#SBATCH -t 48:00:00
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -p batch
#SBATCH --mem=150g
#SBATCH -o GenomeWiseCalculations.out
#SBATCH -e GenomeWiseCalculations.err
#SBATCH --account=epscor-condo

module load samtools/1.16.1-txuglks 
module load r/4.2.2-z6qdiis
module load bcftools/1.16-ewu6fpe
module load htslib/1.17-zxcat2k
module load libnsl

freebayes -f Epsilon_pangenome.fasta -L bam.list -v Epsilon.Freebayes.noMAF.vcf -F 0 --use-best-n-alleles 10 -d -C 1 -p 1 --pooled-continuous -g 1000 -m 30 -q 20 --min-coverage 10 --haplotype-length 0 --report-monomorphic

cat Epsilon.Freebayes.noMAF.vcf | vcf-sort -c > Epsilon.Freebayes.noMAF.sorted.vcf
bgzip -c Epsilon.Freebayes.noMAF.sorted.vcf > Epsilon.Freebayes.noMAF.vcf.gz
tabix -p vcf Epsilon.Freebayes.noMAF.vcf.gz
bcftools filter -g 5 -o Epsilon.Freebayes.minFiltered.vcf -O v Epsilon.Freebayes.noMAF.vcf.gz
vcftools --vcf Epsilon.Freebayes.minFiltered.vcf --max-missing 0.75 --recode --recode-INFO-all --out Epsilon.Freebayes.minFiltered
vcftools --vcf Epsilon.Freebayes.minFiltered.recode.vcf --missing-indv --out Epsilon.Freebayes.minFiltered.recode

cat Epsilon.Freebayes.FINAL.recode.vcf | vcf-sort -c > Epsilon.Freebayes.FINAL.recode.sorted.vcf
bgzip -c Epsilon.Freebayes.FINAL.recode.sorted.vcf > Epsilon.Freebayes.FINAL.recode.sorted.vcf.gz
tabix -p vcf Epsilon.Freebayes.FINAL.recode.sorted.vcf.gz

for file in `cat filelist.txt`
do
bcftools view -e 'REF="N"' -m2 -M2 -v snps -U -s ${file} -Ov -o ${file}.FINAL.reformat.vcf Epsilon.Freebayes.FINAL.recode.sorted.vcf.gz
grep -v "^" ${file}.FINAL.reformat.vcf | perl -anle 'print $F[0] . "\t" . $F[1] if $F[3] =~ /[A-Z]{2}/' > regions.txt
grep -wvf regions.txt ${file}.FINAL.reformat.vcf > tmp
mv tmp ${file}.FINAL.reformat.vcf
done

eval "$(conda shell.bash hook)"
conda activate python2
python /gpfs/data/rbeinart/cbreusing/Scripts/vcf_to_mergedcounts.py .
conda deactivate

for v in Before After ELSC
do
s=`head -n 1 ${v}.list`
awk '{print $1"\t"$2"\t"$3}' ${s}_SNPs.list > ${v}_SNPs.list
  for file in `cat ${v}.list`
  do
  awk '{print $4"\t"$5"\t"$6"\t"$7}' ${file}_SNPs.list | paste ${v}_SNPs.list - > tmp
  mv tmp ${v}_SNPs.list
  done
perl /gpfs/data/rbeinart/cbreusing/Scripts/count_col_entries.pl ${v}_SNPs.list > ${v}_SNPs.merged.list
echo -e "${v}\t${v}_SNPs.merged.list" >> Epsilon.merged.txt
done

conda activate python2
# Calculation for samples pooled by population
Rscript /gpfs/data/rbeinart/cbreusing/Scripts/structure.r Epsilon.merged.txt
python /gpfs/data/rbeinart/cbreusing/Scripts/genome_wise_calculations.py Epsilon.merged.txtFst_pos.txt Epsilon_pangenome.fasta Epsilon.merged

# Calculation for individual samples (to obtain inter-host pi values)
Rscript /gpfs/data/rbeinart/cbreusing/Scripts/structure.r Epsilon.txt
python /gpfs/data/rbeinart/cbreusing/Scripts/genome_wise_calculations.py Epsilon.txtFst_pos.txt Epsilon_pangenome.fasta Epsilon

conda deactivate

vcftools --vcf Epsilon.Freebayes.minFiltered.recode.vcf --haploid --TajimaD 500 --out Before --keep Before.list
vcftools --vcf Epsilon.Freebayes.minFiltered.recode.vcf --haploid --TajimaD 500 --out After --keep After.list
