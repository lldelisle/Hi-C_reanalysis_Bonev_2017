#!/bin/bash

#SBATCH -o slurm-%x-%A_%2a.out # Template for the std output of the job uses the job name, the job id and the array id
#SBATCH -e slurm-%x-%A_%2a.err # Template for the std error of the job
#SBATCH --nodes 1 # We always use 1 node
#SBATCH --ntasks 1 # In this script everything is sequencial
#SBATCH --mem 500G # The memory needed depends on the size of the genome and the size of fastqs
#SBATCH --cpus-per-task 1
#SBATCH --time 72:00:00 # This depends on the size of the fastqs
#SBATCH --array=1-12 # Put here the rows from the table that need to be processed in the table
#SBATCH --job-name HiC_Bonev_dedup # Job name that appear in squeue as well as in output and error text files
#SBATCH --chdir /scratch/ldelisle/HiC_Bonev_mm39/ # This directory must exist, this is where will be the error and out files
## Specific to jed:
#SBATCH --qos serial

# This script run hicup_deduplicator
# convert hicup bam to juicebox format valid pairs
# If needed will subset to pairs both in capture region
# Will sort the pairs to build matrices in cool format
# Balance with cooler
# Do some plots with pyGenomeTracks


##################################
#### TO SET FOR EACH ANALYSIS ####
##################################

### Specify the options for your analysis:
genome=mm39
# For classical DpnII HiC
restName="DpnII"
restOption="--re1 ^GATC,DpnII"
# For Arima protocol
# restName="DpnII_Arima"
# restOption="--arima" # Equivalent to --re1 ^GATC,DpnII:G^ANTC,Arima
# Specify the bin size of matrices to generate (in kb) separated by space
bins="20 5"
# Define a test region for a pgt plot (must be inside the captured region if it is a capture)
# chr7:155000000-158000000 SHH hg38
# chr2:174800000-177800000 HOXD hg38
# chr3:65103500-68603411 Shox2 CaptureC mm39
# chr2:73150000-76150000 HoxD mm39
# chr2:73779626-75669724 HoxD mm10
testRegion="chr2:73779626-75669724"
# Captured region, if it is not a capture, comment these 3 lines
# HoxD mm10:
# chrUsed="chr2"
# start="72402000"
# end="77000000"
# generatefullpairfile="Y" # Comment this line if you only want the file for the cature region

### Specify the paths

# Put in dirPathWithResults the directory
# where the digester file will be stored
# and where a directory will be created
# for each sample
dirPathWithResults="$PWD/"

# All samples are registered into a table where
# first column is the sample name
# second column is the path of the results to deduplicate separated by comma
filePathForTable="/home/ldelisle/softwares/Hi-C_reanalysis_Bonev_2017/scripts_mm39/HiC_Bonev_GSM.txt"
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

dirPathWithCustomHicup2juicer="/home/ldelisle/softwares/Hi-C_reanalysis_Bonev_2017/scripts_mm39/"

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
v=$(hicup --version)
if [ $? -ne 0 ]
then
  echo "HiCUP is not installed but required. Please install it from github (just untar) or conda"
  exit 1
fi
echo $v >> ${fileWithVersions}
v=$(${dirPathWithCustomHicup2juicer}hicup2juicer --version)
if [ $? -ne 0 ]
then
  echo "hicup2juicer is not installed but required. Please install it from github (just untar) or conda"
  exit 1
fi
echo $v >> ${fileWithVersions}
v=$(samtools --version)
if [ $? -ne 0 ]
then
  echo "samtools is not installed but required. Please install it for example in the conda environment."
  exit 1
fi
echo $v >> ${fileWithVersions}
v=$(R --version)
if [ $? -ne 0 ]
then
  echo "R is not installed but required. Please install it for example in the conda environment."
  exit 1
fi
echo $v >> ${fileWithVersions}
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

# Define the output of HiCUP Digester
pathForDigest="${dirPathWithResults}/${genome}_digester_${restName}.txt.gz"
# Define the file with chromosome size (it must correspond to the fasta)
filePathForSizes="${filePathForFasta}.fai"

if [ ! -e $pathForDigest ]; then
  # In order to run the digestion only once
  # This block is only executed for the first sample
  if [ ${SLURM_ARRAY_TASK_ID} = ${SLURM_ARRAY_TASK_MIN} ]; then
    hicup_digester ${restOption} --genome ${genome} --zip --outdir ${dirPathWithResults} ${filePathForFasta}
    mv ${dirPathWithResults}/Digest_${genome}_${restName}* ${pathForDigest}
  else
    echo "Waiting for the first sample to generate the digester file."
    i=0
    while [[ "$i" -lt 20 ]]; do
      sleep 1m
      if [ -e $pathForDigest ]; then
        break
      fi
      ((i++))
    done
    if [ ! -e $pathForDigest ]; then
      echo "After 20 minutes the digester file was not created"
      echo "Checkout what happened or generate it before running"
      exit 1
    fi
  fi
fi

# Get the sample name and intermediate sample names from the table
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

# Only run hicup_deduplicator if the bam does not exists
if [ ! -e merged.dedup.bam ]; then
  if [ ! -e merged.bam ]; then

    # The path of intermediate bams is computed
    intermediate_bams=""
    for s in $(echo $inputs_comma | tr "," " "); do
      intermediate_bams="${intermediate_bams} $(ls ${dirPathWithResults}/${s}/*hicup.bam)"
    done

    echo "BAM to merge are"
    echo $intermediate_bams

    # Merge all intermediate hicup results
    samtools merge -o merged.bam ${intermediate_bams}
  fi
  # Run hicup_deduplicator
  hicup_deduplicator --zip merged.bam
else
  echo "Dedup bam already exists. Not regenerating it."
fi
# if [ ! -e ${sample}.validPairs_nonSorted.txt.gz ]; then
#   hicup2juicer --digest $pathForDigest --zip --useMid merged.dedup.bam
#   mv *prejuicer.gz ${sample}.validPairs_nonSorted.txt.gz
# else
#   echo "Valid pair file already generated. Will not regenerate it."
# fi
if [ ! -e ${sample}.validPairs_nonSorted_new.txt.gz ]; then
  ${dirPathWithCustomHicup2juicer}hicup2juicer --digest $pathForDigest --zip --useMid --nosort merged.dedup.bam
  mv *prejuicer.gz ${sample}.validPairs_nonSorted_new.txt.gz
else
  echo "Valid pair file already generated. Will not regenerate it."
fi

# inputValidPairs=${sample}.validPairs_nonSorted_new.txt.gz

# if [ -z $chrUsed ]; then
#   pairs=${sample}.validPairs.csort.gz
# else
#   if [ -z ${generatefullpairfile} ]; then
#     inputValidPairs_onCapture=${sample}.validPairs_onCapture.txt.gz
#     pairs=${sample}_onCapture.validPairs.csort.gz
#     if [ ! -e ${inputValidPairs_onCapture} ]; then
#         zcat ${inputValidPairs} | awk -v c=${chrUsed} -v s=${start} -v e=${end} '$3==c && $7==c && $4>s && $4<e && $8>s && $8<e{print}' | gzip > ${inputValidPairs_onCapture}
#     else
#       echo "Valid pairs on capture already generated. Will not regenerate it."
#     fi
#     # Update the inputValidPairs variable:
#     inputValidPairs=${inputValidPairs_onCapture}
#   else
#     pairs=${sample}.validPairs.csort.gz
#   fi
# fi

# if [ ! -e $filePathForSizes ]; then
#   samtools faidx $filePathForFasta
# fi

# if [ ! -e $pairs ]; then
#   # sort and index the pairs with cooler and tabix
#   cooler csort -i tabix -c1 3 -c2 7 -p1 4 -p2 8 -o $pairs $inputValidPairs $filePathForSizes
# fi

# if [ ! -z $chrUsed ] && [ ! -z ${generatefullpairfile} ]; then
#   # Extract capture region from global sorted pair:
#   tabix $pairs ${chrUsed}:${start}-${end} | awk -v c=${chrUsed} -v s=${start} -v e=${end} '$7==c && $8>s && $8<e{print}' | bgzip > ${sample}_onCapture.validPairs.csort.gz
#   pairs=${sample}_onCapture.validPairs.csort.gz
#   tabix -s 3 -b 4 -e 4 ${pairs}
# fi

# ini_file=${sample}.ini
# echo "[x-axis]" > $ini_file

# for bin in $bins; do
#   if [ ! -e ${genome}.${bin}kb.bins ]; then
#     cooler makebins $filePathForSizesForBin "${bin}000" > ${genome}.${bin}kb.bins
#   fi
#   if [ ! -e ${sample}.${bin}kb.cool ]; then
#     echo "cload $pairs"
#     # This failed see https://github.com/open2c/cooler/issues/315
#     # cooler cload tabix -p 1 -c2 7 -p2 8 --assembly $genome ${genome}.${bin}kb.bins $pairs ${sample}.${bin}kb.cool
#     cooler cload tabix -c2 7 -p2 8 --assembly $genome ${genome}.${bin}kb.bins $pairs ${sample}.${bin}kb.cool
#     cp ${sample}.${bin}kb.cool ${sample}_raw.${bin}kb.cool
#     echo "Balancing"
#     cooler balance --cis-only ${sample}.${bin}kb.cool
#     echo "Balanced"
#   fi
#   echo "[${sample}_${bin}kb]
# file = ${sample}.${bin}kb.cool
# depth = 2000000
# min_value = 0
# title = ${sample}_${bin}kb
# [spacer]
# [${sample}_raw_${bin}kb]
# file = ${sample}_raw.${bin}kb.cool
# depth = 2000000
# min_value = 0
# title = ${sample}_raw_${bin}kb
# [spacer]
# " >> ${ini_file}
# done
# echo "[x-axis]" >> ${ini_file}
# # Generate a basic plot on the testRegion
# pgt --tracks ${ini_file} --region ${testRegion} --fontSize 6 -o ${ini_file/.ini/_testRegion.pdf}

# if [ ! -z $chrUsed ]; then
#   # If it is a capture, plot on the whole capture region:
#   pgt --tracks ${ini_file} --region ${chrUsed}:${start}-${end} --fontSize 6 -o ${ini_file/.ini/.pdf}
# fi

# # Copy all final files to specific directories
# mkdir -p ${dirPathWithResults}/allFinalFiles/cool
# cp *.cool ${dirPathWithResults}/allFinalFiles/cool/
# mkdir -p ${dirPathWithResults}/allFinalFiles/pairs
# cp ${pairs} ${dirPathWithResults}/allFinalFiles/pairs/
# mkdir -p ${dirPathWithResults}/allFinalFiles/visualisation
# cp *.pdf ${dirPathWithResults}/allFinalFiles/visualisation

