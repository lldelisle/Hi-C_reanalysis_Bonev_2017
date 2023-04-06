#!/bin/bash

#SBATCH -o slurm-%x-%A_%2a.out # Template for the std output of the job uses the job name, the job id and the array id
#SBATCH -e slurm-%x-%A_%2a.err # Template for the std error of the job
#SBATCH --nodes 1 # We always use 1 node
#SBATCH --ntasks 1 # In this script everything is sequencial
#SBATCH --mem 150G # The memory needed depends on the size of the genome and the size of fastqs
#SBATCH --cpus-per-task 1
#SBATCH --time 72:00:00 # This depends on the size of the fastqs
#SBATCH --array=1-3 # Put here the rows from the table that need to be processed in the table
#SBATCH --job-name HiC_Bonev_merge # Job name that appear in squeue as well as in output and error text files
#SBATCH --chdir /scratch/ldelisle/HiC_Bonev_mm39/ # This directory must exist, this is where will be the error and out files
## Specific to jed:
#SBATCH --qos serial

# This script will combine all pairs from replicates
# Will sort the pairs per chromosome and build matrices in cool format for chr2
# Balance with cooler
# Do some plots with pyGenomeTracks


##################################
#### TO SET FOR EACH ANALYSIS ####
##################################

### Specify the options for your analysis:
genome=mm39
# Specify the bin size of matrices to generate (in kb) separated by space
bins="5"
# Define a test region for a pgt plot (must be inside the captured region if it is a capture)
# chr7:155000000-158000000 SHH hg38
# chr2:174800000-177800000 HOXD hg38
# chr3:65103500-68603411 Shox2 CaptureC mm39
# chr2:73150000-76150000 HoxD mm39
# chr2:73779626-75669724 HoxD mm10
testRegion="chr2:73779626-75669724"

### Specify the paths

# Put in dirPathWithResults the directory
# where the digester file will be stored
# and where a directory will be created
# for each sample
dirPathWithResults="$PWD/"

# All samples are registered into a table where
# first column is the sample name
# second column is the path of the results to deduplicate separated by comma
filePathForTable="/home/ldelisle/softwares/Hi-C_reanalysis_Bonev_2017/scripts_mm39/HiC_Bonev_samples.txt"
## The script will write a file with software versions close to this file
fileWithVersions=$(dirname $filePathForTable)/${SLURM_JOBID}_versions.txt

filePathForFasta="/home/ldelisle/genomes/fasta/${genome}.fa"
# If you downloaded HiCUP from github
# Put here the directory it will be added to the PATH:
# dirPathForHiCUP="/home/ldelisle/softwares/HiCUP-0.9.2/"
# If you used conda, just comment it

# You can decide to generate cool files only for part of the genome
# Use it if you want to bin only specific chromosome (ie exclude contigs)
# filePathForSizesForBin="/home/users/d/darbellf/live/genomes/${genome}/size_my_favorite_chr.txt"
filePathForSizesForBin="${filePathForFasta}.fai"

### Specify the way to deal with dependencies:

# conda create -n hic202302 python=3.10 mamba
# conda activate hic202302
# mamba install -c bioconda -c conda-forge pygenometracks 'hicup>=0.9.2'
condaEnvName=hic202302


##################################
####### BEGINING OF SCRIPT #######
##################################

if [ ! -z ${dirPathForHiCUP} ]; then
  export PATH=$PATH:${dirPathForHiCUP}:${dirPathForHiCUP}/Conversion/
fi
# Check everything is set correctly:
# Conda environment:
# This line is to adapt the conda to the shell
source $(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh
# We check if the conda environment exists
exists=$(conda info --envs | awk -v ce=${condaEnvName} '$1==ce{print}' | wc -l)
# It if does not exists an error is raised
if [ $exists -ne 1 ]; then
  echo "conda environment ${condaEnvName} does not exists. Create it before."
  exit 1
fi
# Activate the conda environment
conda activate ${condaEnvName}
# Check all softwares are present and write version to stdout:
v=$(cooler --version)
# Check cooler is installed:
if [ $? -ne 0 ]
then
  echo "Cooler is not installed but required. Please install it for example in the conda environment"
  exit 1
fi
echo $v >> ${fileWithVersions}
v=$(pgt --version)
# Check pgt is installed:
if [ $? -ne 0 ]
then
  echo "pyGenomeTracks is not installed but required. Please install it for example in the conda environment"
  exit 1
fi
echo $v >> ${fileWithVersions}

# Define the file with chromosome size (it must correspond to the fasta)
filePathForSizes="${filePathForFasta}.fai"

# Get the sample name and replicate sample names from the table
sample=$(cat ${filePathForTable} | awk -v i=${SLURM_ARRAY_TASK_ID} 'NR==i{print $1}')
inputs_comma=$(cat ${filePathForTable} | awk -v i=${SLURM_ARRAY_TASK_ID} 'NR==i{print $2}')

# Each sample is processed into an independent directory:
pathResults=${dirPathWithResults}/${sample}/

# The directory is created (if not existing)
mkdir -p ${pathResults}

# The name of the sample is written in stdout
echo ${sample}

# The analysis part takes part within the pathResults
cd ${pathResults}

# Concatenate the non sorted pairs
if [ ! -e ${sample}.validPairs_notSorted.txt.gz ]; then
  # The path of intermediate validPairs is computed
  intermediate_pairs=""
  for s in $(echo $inputs_comma | tr "," " "); do
      # intermediate_pairs="${intermediate_pairs} ${s}.validPairs_nonSorted.txt.gz"
      intermediate_pairs="${intermediate_pairs} ${dirPathWithResults}/${s}/${s}.validPairs_nonSorted_new.txt.gz"
  done

  echo "Pairs to merge are"
  echo $intermediate_pairs

  # Merge all intermediate pairs
  cat ${intermediate_pairs} > ${sample}.validPairs_notSorted.txt.gz
else
  echo "Valid pair file already generated. Will not regenerate it."
fi

# Split by chromosome for cis and put all trans together
nbChrGZ=$(ls ${sample}_split_chr*validPairs_notSorted.txt.gz | wc -l)

if [ $nbChrGZ -eq 0 ]; then
  zcat ${sample}.validPairs_notSorted.txt.gz | awk -v base="${sample}_split" '$3==$7{print > base"_"$3".validPairs_notSorted.txt"}$3!=$7{print}' | gzip > ${sample}_split_trans.validPairs_notSorted.txt.gz
  gzip ${sample}_split_chr*.validPairs_notSorted.txt
fi

if [ ! -e $filePathForSizes ]; then
  samtools faidx $filePathForFasta
fi


##NEEDED FOR BIG FILES##
mkdir -p temp
export TMPDIR=$PWD/temp


mkdir -p logs

for f in ${sample}_split_*.validPairs_notSorted.txt.gz; do
  tname=$(basename $f .validPairs_notSorted.txt.gz)
  tpairs=${tname}.validPairs.csort.gz
  if [ ! -e ${tpairs}.tbi ]; then
    cooler csort -i tabix -c1 3 -c2 7 -p1 4 -p2 8 -o $tpairs $f $filePathForSizes 2> logs/${tname}.log
  fi
done

ini_file=${sample}.ini
echo "[x-axis]" > $ini_file

for bin in $bins; do
  if [ ! -e ${genome}.${bin}kb.bins ]; then
    cooler makebins $filePathForSizesForBin "${bin}000" > ${genome}.${bin}kb.bins
  fi
  if [ ! -e ${sample}.${bin}kb.cool ]; then
    cooler cload tabix -c2 7 -p2 8 --assembly $genome ${genome}.${bin}kb.bins ${sample}_split_chr2.validPairs.csort.gz ${sample}_chr2.${bin}kb.cool
    cp ${sample}_chr2.${bin}kb.cool ${sample}_chr2_raw.${bin}kb.cool
    echo "Balancing"
    cooler balance --cis-only ${sample}_chr2.${bin}kb.cool
    echo "Balanced"
  fi
  echo "[${sample}_${bin}kb]
file = ${sample}_chr2.${bin}kb.cool
depth = 2000000
min_value = 0
title = ${sample}_${bin}kb
[spacer]
[${sample}_raw_${bin}kb]
file = ${sample}_chr2_raw.${bin}kb.cool
depth = 2000000
min_value = 0
title = ${sample}_raw_${bin}kb
[spacer]
" >> ${ini_file}
done
echo "[x-axis]" >> ${ini_file}
# Generate a basic plot on the testRegion
pgt --tracks ${ini_file} --region ${testRegion} --fontSize 6 -o ${ini_file/.ini/_testRegion.pdf}


# Copy all final files to specific directories
mkdir -p ${dirPathWithResults}/allFinalFiles/cool
cp *.cool ${dirPathWithResults}/allFinalFiles/cool/
mkdir -p ${dirPathWithResults}/allFinalFiles/pairs
cp *csort.gz ${dirPathWithResults}/allFinalFiles/pairs/
mkdir -p ${dirPathWithResults}/allFinalFiles/visualisation
cp *.pdf ${dirPathWithResults}/allFinalFiles/visualisation

# For GEO
mkdir ${dirPathWithResults}/toGEO
cp *csort.gz ${dirPathWithResults}/toGEO/
