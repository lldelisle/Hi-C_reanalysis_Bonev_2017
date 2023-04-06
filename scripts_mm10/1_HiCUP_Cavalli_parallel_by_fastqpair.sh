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
#SBATCH --array=1-102
#SBATCH --job-name HiCUP
#SBATCH --workdir /scratch/ldelisle/HiCUP

path="$PWD/"
pathForScripts="/home/ldelisle/softwares/Hi-C_reanalysis_Bonev_2017/scripts/"
pathForTableWithSRA="/home/ldelisle/softwares/Hi-C_reanalysis_Bonev_2017/table.txt"
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


pathForDigest="${path}/${genome}_digester_$restName.txt.gz"

#Digest 
if [ $SLURM_ARRAY_TASK_ID = 1 ]; then
  if [ ! -e $pathForDigest ]; then
    ${pathForHiCUP}hicup_digester --re1 $restSeq,$restName --genome $genome --zip --outdir $path $pathForFasta
    mv ${path}/Digest_${genome}_${restName}* $pathForDigest
  fi
fi


sample=`cat $pathForTableWithSRA | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $2}'`

mkdir -p fromGEO
module load sra-toolkit/2.9.2
bash ${pathForScripts}/getPEFastqFromSRAtableAndIndex_fast_parallel_fasterq_secured.sh $pathForTableWithSRA $SLURM_ARRAY_TASK_ID

fastqFileR1="${path}/fromGEO/${sample}_R1.fastq.gz"
fastqFileR2="${path}/fromGEO/${sample}_R2.fastq.gz"

mkdir -p $sample
cd $sample

${pathForHiCUP}hicup --bowtie2 $pathForBowtie2 --digest $pathForDigest --format Sanger --index $pathForB2Index --keep --threads 28 --zip --r $pathForR $fastqFileR1 $fastqFileR2
