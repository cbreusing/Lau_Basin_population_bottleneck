!/bin/bash
SBATCH -J SNPCalling
SBATCH -t 200:00:00
SBATCH -n 8
SBATCH -N 1
SBATCH -p batch
SBATCH --mem=50g
SBATCH -o SNPCalling-BW.out
SBATCH -e SNPCalling-BW.err
SBATCH --account=epscor-condo

module load bowtie2/2.4.2-xdquyzq
module load samtools/1.16.1-txuglks
module load gatk/4.3.0.0-234wqft 
module load picard/2.26.2-qabtyqy
module load r/4.2.2-z6qdiis
module load fastme/2.1.5.1-kmg5til 
module load gsl/2.7.1-khmyfcy
module load libnsl/2.0.1-ed2i5hn
module load angsd/0.935-cbhuwc7

for file in `cat filelist.txt`
do
bowtie2 --very-sensitive-local -p ${SLURM_CPUS_ON_NODE} -x A_boucheti.Trinity.merged95.filtered.fasta -1 ../../${file}_R1_clean.fastq -2 ../../${file}_R2_clean.fastq | samtools view -bS -h -@ ${SLURM_CPUS_ON_NODE} - | samtools sort -@ ${SLURM_CPUS_ON_NODE} - > ${file}.sorted.bam
picard MarkDuplicates I=${file}.sorted.bam O=${file}.dedup.bam M=${file}.metrics.txt
samtools index ${file}.dedup.bam
eval "$(conda shell.bash hook)"
conda activate lofreq
lofreq viterbi -f A_boucheti.Trinity.merged95.filtered.fasta -k ${file}.dedup.bam | samtools sort -@ ${SLURM_CPUS_ON_NODE} - > ${file}.realigned.bam
lofreq indelqual -f A_boucheti.Trinity.merged95.filtered.fasta --dindel -o ${file}.indelqual.bam ${file}.realigned.bam
conda deactivate
samtools index ${file}.indelqual.bam
done

rm *.sorted.bam
rm *.dedup.bam*
rm *.realigned.bam

ls *indelqual.bam > bam.list

angsd -P ${SLURM_CPUS_ON_NODE} -bam bam.list -ref A_boucheti.Trinity.merged95.filtered.fasta -out ANGSD_basin-wide/A_boucheti.qc -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -C 50 -baq 1 -minMapQ 30 -doQsDist 1 -doDepth 1 -doCounts 1 -maxDepth 1000
Rscript /gpfs/data/rbeinart/Software/ngsTools/Scripts/plotQC.R ANGSD_basin-wide/A_boucheti.qc

for POP in KM TC TM ABE;
do
IND=`cat ANGSD_basin-wide/${POP}.list | wc -l`
FRAC=`echo "${IND}*0.5" | bc | awk '{print int($1+0.5)}'`
angsd -P ${SLURM_CPUS_ON_NODE} -bam ANGSD_basin-wide/${POP}.list -ref A_boucheti.Trinity.merged95.filtered.fasta -gl 1 -baq 1 -C 50 -minInd ${FRAC} -minMapQ 30 -minQ 20 -setMinDepth 2 -setMaxDepth 160 -doCounts 1 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -doMaf 1 -doMajorMinor 1 -SNP_pval 1e-6 -doGlf 3 -minMaf 0.01 -skipTriallelic 1 -out ANGSD_basin-wide/${POP}
zcat ANGSD_basin-wide/${POP}.glf.pos.gz > ANGSD_basin-wide/${POP}.glf.pos
cut -f1 ANGSD_basin-wide/${POP}.glf.pos | sort | uniq > ANGSD_basin-wide/${POP}.chrs.txt
grep -f ANGSD_basin-wide/${POP}.chrs.txt ANGSD_basin-wide/${POP}.glf.pos | awk -F"\t" '!_[$1]++' | perl -anle 'print $F[0] . "\t" . $F[1]' > ANGSD_basin-wide/${POP}.unlinked.pos
angsd sites index ANGSD_basin-wide/${POP}.unlinked.pos
angsd -P ${SLURM_CPUS_ON_NODE} -bam ANGSD_basin-wide/${POP}.list -ref A_boucheti.Trinity.merged95.filtered.fasta -gl 1 -baq 1 -C 50 -minInd ${FRAC} -sites ANGSD_basin-wide/${POP}.unlinked.pos -rf ANGSD_basin-wide/${POP}.chrs.txt -minMapQ 30 -minQ 20 -setMinDepth 2 -setMaxDepth 160 -doCounts 1 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -doMaf 1 -doMajorMinor 1 -SNP_pval 1e-6 -doGlf 3 -minMaf 0.01 -skipTriallelic 1 -out ANGSD_basin-wide/${POP}
NSITES=`zcat ANGSD_basin-wide/${POP}.mafs.gz | tail -n+2 | wc -l`
zcat ANGSD_basin-wide/${POP}.glf.gz > ANGSD_basin-wide/${POP}.glf
/gpfs/data/rbeinart/Software/ngsTools/ngsF/ngsF.sh --n_ind ${IND} --n_sites ${NSITES} --glf ANGSD_basin-wide/${POP}.glf --out ANGSD_basin-wide/${POP}.indF
done

head -n 10 ANGSD_basin-wide/ABE.indF > ANGSD_basin-wide/ABE_B.indF
tail -n 20 ANGSD_basin-wide/ABE.indF > ANGSD_basin-wide/ABE_A.indF
cat ANGSD_basin-wide/KM.indF ANGSD_basin-wide/TC.indF ANGSD_basin-wide/ABE_B.indF > ANGSD_basin-wide/Before.indF
cat ANGSD_basin-wide/TM.indF ANGSD_basin-wide/ABE_A.indF > ANGSD_basin-wide/After.indF
cat ANGSD_basin-wide/After.indF ANGSD_basin-wide/Before.indF > ANGSD_basin-wide/ELSC.indF

POP=ELSC
IND=`cat ${POP}.list | wc -l`
FRAC=`echo "${IND}*0.75" | bc | awk '{print int($1+0.5)}'`
angsd -P ${SLURM_CPUS_ON_NODE} -bam ${POP}.list -ref A_boucheti.Trinity.merged95.filtered.fasta -gl 1 -doGlf 1 -baq 1 -C 50 -minInd ${FRAC} -minMapQ 30 -minQ 20 -setMinDepth 2 -setMaxDepth 160 -doCounts 1 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -doMaf 1 -doMajorMinor 1 -SNP_pval 1e-6 -dosnpstat 1 -doHWE 1 -sb_pval 0.05 -hetbias_pval 0.05 -makeMatrix 1 -doIBS 1 -doCov 1 -doPost 1 -doGeno 2 -minMaf 0.01 -indF ANGSD_basin-wide/${POP}.indF -skipTriallelic 1 -out ANGSD_basin-wide/A_boucheti
zcat ANGSD_basin-wide/${POP}.mafs.gz | tail -n+2 | perl -anle 'print $F[0] . "\t" . $F[1]' > ANGSD_basin-wide/sites.txt

for POP in KM TC TM ABE;
do
IND=`cat ANGSD_basin-wide/${POP}.list | wc -l`
FRAC=`echo "${IND}*0.5" | bc | awk '{print int($1+0.5)}'`
angsd -P ${SLURM_CPUS_ON_NODE} -bam ANGSD_basin-wide/${POP}.list -ref A_boucheti.Trinity.merged95.filtered.fasta -anc A_boucheti.Trinity.merged95.filtered.fasta -gl 1 -baq 1 -C 50 -minInd ${FRAC} -minMapQ 30 -minQ 20 -setMinDepth 2 -setMaxDepth 160 -doCounts 1 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -doMaf 1 -doMajorMinor 4 -doSaf 2 -indF ANGSD_basin-wide/${POP}.indF -out ANGSD_basin-wide/${POP}
realSFS ANGSD_basin-wide/${POP}.saf.idx -tole 1e-6 -maxIter 5000 -P ${SLURM_CPUS_ON_NODE} -fold 1 > ANGSD_basin-wide/${POP}.sfs
done

realSFS ANGSD_basin-wide/KM.saf.idx ANGSD_basin-wide/TC.saf.idx -tole 1e-6 -maxIter 5000 -fold 1 -P ${SLURM_CPUS_ON_NODE} > ANGSD_basin-wide/KM-TC.folded.sfs
realSFS fst index ANGSD_basin-wide/KM.saf.idx ANGSD_basin-wide/TC.saf.idx -sfs ANGSD_basin-wide/KM-TC.folded.sfs -fold 1 -fstout ANGSD_basin-wide/KM-TC -whichFST 1 -maxIter 5000
realSFS ANGSD_basin-wide/KM.saf.idx ANGSD_basin-wide/ABE.saf.idx -tole 1e-6 -maxIter 5000 -fold 1 -P ${SLURM_CPUS_ON_NODE} > ANGSD_basin-wide/KM-ABE.folded.sfs
realSFS fst index ANGSD_basin-wide/KM.saf.idx ANGSD_basin-wide/ABE.saf.idx -sfs ANGSD_basin-wide/KM-ABE.folded.sfs -fold 1 -fstout ANGSD_basin-wide/KM-ABE -whichFST 1 -maxIter 5000
realSFS ANGSD_basin-wide/KM.saf.idx ANGSD_basin-wide/TM.saf.idx -tole 1e-6 -maxIter 5000 -fold 1 -P ${SLURM_CPUS_ON_NODE} > ANGSD_basin-wide/KM-TM.folded.sfs
realSFS fst index ANGSD_basin-wide/KM.saf.idx ANGSD_basin-wide/TM.saf.idx -sfs ANGSD_basin-wide/KM-TM.folded.sfs -fold 1 -fstout ANGSD_basin-wide/KM-TM -whichFST 1 -maxIter 5000
realSFS ANGSD_basin-wide/TC.saf.idx ANGSD_basin-wide/ABE.saf.idx -tole 1e-6 -maxIter 5000 -fold 1 -P ${SLURM_CPUS_ON_NODE} > ANGSD_basin-wide/TC-ABE.folded.sfs
realSFS fst index ANGSD_basin-wide/TC.saf.idx ANGSD_basin-wide/ABE.saf.idx -sfs ANGSD_basin-wide/TC-ABE.folded.sfs -fold 1 -fstout ANGSD_basin-wide/TC-ABE -whichFST 1 -maxIter 5000
realSFS ANGSD_basin-wide/TC.saf.idx ANGSD_basin-wide/TM.saf.idx -tole 1e-6 -maxIter 5000 -fold 1 -P ${SLURM_CPUS_ON_NODE} > ANGSD_basin-wide/TC-TM.folded.sfs
realSFS fst index ANGSD_basin-wide/TC.saf.idx ANGSD_basin-wide/TM.saf.idx -sfs ANGSD_basin-wide/TC-TM.folded.sfs -fold 1 -fstout ANGSD_basin-wide/TC-TM -whichFST 1 -maxIter 5000
realSFS ANGSD_basin-wide/ABE.saf.idx ANGSD_basin-wide/TM.saf.idx -tole 1e-6 -maxIter 5000 -fold 1 -P ${SLURM_CPUS_ON_NODE} > ANGSD_basin-wide/ABE-TM.folded.sfs
realSFS fst index ANGSD_basin-wide/ABE.saf.idx ANGSD_basin-wide/TM.saf.idx -sfs ANGSD_basin-wide/ABE-TM.folded.sfs -fold 1 -fstout ANGSD_basin-wide/ABE-TM -whichFST 1 -maxIter 5000

realSFS fst stats ANGSD_basin-wide/KM-TC.fst.idx -whichFST 1 -maxIter 5000 > ANGSD_basin-wide/KM-TC.fst.txt
realSFS fst stats ANGSD_basin-wide/KM-ABE.fst.idx -whichFST 1 -maxIter 5000 > ANGSD_basin-wide/KM-ABE.fst.txt
realSFS fst stats ANGSD_basin-wide/KM-TM.fst.idx -whichFST 1 -maxIter 5000 > ANGSD_basin-wide/KM-TM.fst.txt
realSFS fst stats ANGSD_basin-wide/TC-ABE.fst.idx -whichFST 1 -maxIter 5000 > ANGSD_basin-wide/TC-ABE.fst.txt
realSFS fst stats ANGSD_basin-wide/TC-TM.fst.idx -whichFST 1 -maxIter 5000 > ANGSD_basin-wide/TC-TM.fst.txt
realSFS fst stats ANGSD_basin-wide/ABE-TM.fst.idx -whichFST 1 -maxIter 5000 > ANGSD_basin-wide/ABE-TM.fst.txt

for i in `cat ANGSD_basin-wide/*fst.idx`
do 
realSFS fst print ANGSD_basin-wide/${i} > ANGSD_basin-wide/${i}.out
grep -f ANGSD_basin-wide/sites.txt ANGSD_basin-wide/${i}.out > ANGSD_basin-wide/${i}.sites.out
sum1=`awk -F "\t" '{sum=sum+$3} END{print sum}' ANGSD_basin-wide/${i}.sites.out`
sum2=`awk -F "\t" '{sum=sum+$4} END{print sum}' ANGSD_basin-wide/${i}.sites.out`
echo "scale=4 ; $sum1 / $sum2" | bc > ANGSD_basin-wide/${i}.sites.txt
done

for POP in Before After;
do
IND=`cat ${POP}.list | wc -l`
FRAC=`echo "${IND}*0.75" | bc | awk '{print int($1+0.5)}'`
angsd -P ${SLURM_CPUS_ON_NODE} -bam ${POP}.list -ref A_boucheti.Trinity.merged95.filtered.fasta -anc A_boucheti.Trinity.merged95.filtered.fasta -gl 1 -baq 1 -C 50 -minInd ${FRAC} -minMapQ 30 -minQ 20 -setMinDepth 2 -setMaxDepth 160 -doCounts 1 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -doMaf 1 -doMajorMinor 4 -doSaf 2 -indF ANGSD_basin-wide/${POP}.indF -out ANGSD_basin-wide/${POP}
realSFS ANGSD_basin-wide/${POP}.saf.idx -tole 1e-6 -maxIter 5000 -P ${SLURM_CPUS_ON_NODE} -fold 1 > ANGSD_basin-wide/${POP}.sfs
realSFS saf2theta ANGSD_basin-wide/${POP}.saf.idx -outname ANGSD_basin-wide/${POP} -sfs ANGSD_basin-wide/${POP}.sfs -fold 1
thetaStat do_stat ANGSD_basin-wide/${POP}.thetas.idx
thetaStat print ANGSD_basin-wide/${POP}.thetas.idx > ANGSD_basin-wide/${POP}.thetas
thetaStat do_stat ANGSD_basin-wide/${POP}.thetas.idx -win 500 -step 100 -outnames ANGSD_basin-wide/${POP}.thetas
perl extract_thetas.pl ANGSD_basin-wide/sites.txt ANGSD_basin-wide/${POP}.thetas ANGSD_basin-wide/${POP}.sites.thetas
angsd -P ${SLURM_CPUS_ON_NODE} -bam ${POP}.list -ref A_boucheti.Trinity.merged95.filtered.fasta -gl 1 -baq 1 -C 50 -minInd ${FRAC} -minMapQ 30 -minQ 20 -setMinDepth 2 -setMaxDepth 160 -doCounts 1 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -doMaf 1 -doMajorMinor 1 -doHWE 1 -doGlf 1 -indF ANGSD_basin-wide/${POP}.indF -out ANGSD_basin-wide/${POP}
gunzip ANGSD_basin-wide/${POP}.hwe.gz
echo "Hexp" > ANGSD_basin-wide/${POP}.Hexp
perl -anle 'print 2*($F[4])*(1-$F[4]) if $. > 1' ANGSD_basin-wide/${POP}.hwe >> ANGSD_basin-wide/${POP}.Hexp
paste ANGSD_basin-wide/${POP}.hwe ANGSD_basin-wide/${POP}.Hexp > ANGSD_basin-wide/tmp
echo "Hobs" > ANGSD_basin-wide/${POP}.Hobs
perl -anle 'print $F[9]-$F[6]*$F[9] if $. > 1' ANGSD_basin-wide/tmp >> ANGSD_basin-wide/${POP}.Hobs 
paste ANGSD_basin-wide/tmp ANGSD_basin-wide/${POP}.Hobs > ANGSD_basin-wide/${POP}.hwe
rm ANGSD_basin-wide/tmp
perl extract_hwe.pl ANGSD_basin-wide/sites.txt ANGSD_basin-wide/${POP}.hwe ANGSD_basin-wide/${POP}.sites.hwe
done
