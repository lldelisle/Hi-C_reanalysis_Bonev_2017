#!/bin/bash

#SBATCH -o slurm-%x-%A_%2a.out
#SBATCH -e slurm-%x-%A_%2a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lucille.delisle@epfl.ch
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --mem 20G 
#SBATCH --cpus-per-task 28
#SBATCH --time 72:00:00
#SBATCH --array=1,12,19,26,33,44,51,58,65,68,82,89
#SBATCH --job-name HiCUP_merge
#SBATCH --workdir /scratch/ldelisle/HiCUP

path="$PWD/"
pathForScripts="/home/ldelisle/softwares/Hi-C_reanalysis_Bonev_2017/scripts/"
pathForTableWithSRA="table.txt"
genome="mm10"
restSeq="^GATC"
restName="DpnII"
pathForB2Index="/home/ldelisle/genomes/bowtie2/$genome"
pathForFasta="/home/ldelisle/genomes/fasta/${genome}.fa"
pathForHiCUP="/home/ldelisle/softwares/hicup_v0.6.1/"


module purge
module load intel/18.0.2
module load intel-mkl/2018.2.199
module load r/3.5.0
pathForR=`which R`
module load gcc/6.4.0
module load bowtie2/2.3.4.1
module load samtools/1.8
pathForBowtie2=`which bowtie2`
#For the moment it is not possible to have R and bowtie2

sampleOri=`cat $pathForTableWithSRA | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $2}'`
sample=`basename $sampleOri _1`
mkdir -p $sample

cd $sample

if [ ! -e ${sample}_R1_2.filt.bam ]; then
  # Merge the filter bam
  # Put the header
  samtools view -H ../${sampleOri}/${sampleOri}_R1_2.filt.bam | gzip > ${sample}.filt.sam.gz
  # Then concatenate without the header but gzip to save space
  for fil in ../${sample}_*/*filt.bam; do
    samtools view $fil  | gzip >> ${sample}.filt.sam.gz
  done
  # Convert to bam
  zcat ${sample}.filt.sam.gz | samtools view -b > ${sample}_R1_2.filt.bam
  rm ${sample}.filt.sam.gz
fi

# Remove duplicates
if [ ! -e ${sample}_R1_2.dedup.bam ]; then
  perl ${pathForHiCUP}hicup_deduplicator --r $pathForR --threads 28 --zip ${sample}_R1_2.filt.bam
fi

inputBAM=${sample}_R1_2.dedup.bam

pathForDigestNGZ="${path}/${genome}_digester_$restName.txt"

if [ ! -e $pathForDigestNGZ ]; then
  gunzip -c $pathForDigest > $pathForDigestNGZ
fi

if [ ! -e ${sample}.validPairs_nonSorted.txt.gz ]; then
  python ${pathForScripts}/fromHicupToJuicebox_withTempFolder.py  --fragmentFile $pathForDigestNGZ --colForChr 1 --colForStart 2 --colForEnd 3 --colForID 4 --lineToSkipInFragmentFile 2 --useMid $inputBAM | gzip > ${sample}.validPairs_nonSorted.txt.gz
fi

