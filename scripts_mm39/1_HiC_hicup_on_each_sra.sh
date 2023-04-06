#!/bin/bash

#SBATCH -o slurm-%x-%A_%2a.out # Template for the std output of the job uses the job name, the job id and the array id
#SBATCH -e slurm-%x-%A_%2a.err # Template for the std error of the job
#SBATCH --nodes 1 # We always use 1 node
#SBATCH --ntasks 1 # In this script everything is sequencial
#SBATCH --mem 100G # The memory needed depends on the size of the genome and the size of fastqs
#SBATCH --cpus-per-task 24 # This allows to speed the mapping part of HiCUP
#SBATCH --time 72:00:00 # This depends on the size of the fastqs
#SBATCH --array=1-101 # Put here the rows from the table that need to be processed in the table
#SBATCH --job-name HiC_Bonev # Job name that appear in squeue as well as in output and error text files
#SBATCH --chdir /scratch/ldelisle/HiC_Bonev_mm39/ # This directory must exist, this is where will be the error and out files
## Specific to jed:
#SBATCH --qos serial

# This script run hicup >=0.9.2 with Bowtie2


##################################
#### TO SET FOR EACH ANALYSIS ####
##################################

### Specify the options for your analysis:
# number of CPU to use
# Only change if you don't want to use all CPUs allocated
nbOfThreads=${SLURM_CPUS_PER_TASK}
# Which genome to map on
# /!\ This script will use bowtie2 and bowtie2 is not 'alt-aware' so do not use a genome with alt contigs
genome=mm39
# For classical DpnII HiC
restName="DpnII"
restOption="--re1 ^GATC,DpnII"
# For Arima protocol
# restName="DpnII_Arima"
# restOption="--arima" # Equivalent to --re1 ^GATC,DpnII:G^ANTC,Arima
optionForHiCUP="" # Use "--nofill" if there was no biotin fill-in (ie. CaptureC-like)

### Specify the paths

# Put in dirPathWithResults the directory
# where the digester file will be stored
# and where a directory will be created
# for each sample
dirPathWithResults="$PWD/"
# Where fastqs are stored:
dirPathForFastq="${dirPathWithResults}/fastq/"
# All samples are registered into a table where
# first column is the sample name
# second column is the path of the R1 fastq relatively to dirPathForFastq
# third column is same for R2
# Alternatively second column can be SRA number
filePathForTable="/home/ldelisle/softwares/Hi-C_reanalysis_Bonev_2017/scripts_mm39/HiC_Bonev.txt"
## The script will write a file with software versions close to this file
fileWithVersions=$(dirname $filePathForTable)/${SLURM_JOBID}_versions.txt

basenamePathForB2Index="/scratch/ldelisle/genomes/bowtie2/${genome}"
filePathForFasta="/scratch/ldelisle/genomes/fasta/${genome}.fa"
# If you downloaded HiCUP from github
# Put here the directory it will be added to the PATH:
# dirPathForHiCUP="/home/ldelisle/softwares/HiCUP-0.9.2/"
# If you used conda, just comment it

# conda create -n hic202302 python=3.10 mamba
# conda activate hic202302
# mamba install -c bioconda -c conda-forge pygenometracks 'hicup>=0.9.2'
# If you want to use sra you also need sra-tools>=2.11
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
v=$(hicup --version)
if [ $? -ne 0 ]
then
  echo "HiCUP is not installed but required. Please install it from github (just untar) or conda"
  exit 1
fi
echo $v >> ${fileWithVersions}
v=$(bowtie2 --version)
if [ $? -ne 0 ]
then
  echo "Bowtie2 is not installed but required. Please install it for example in the conda environment."
  exit 1
fi
echo $v >> ${fileWithVersions}
if [ ! -e ${basenamePathForB2Index}.rev.2.bt2 ]; then
  echo "${basenamePathForB2Index}.rev.2.bt2 does not exists. Either your variable basenamePathForB2Index is wrong or your bowtie2 index did not full run."
  exit 1
fi
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
# Check plotly and rmarkdown and tidyverse are installed:
Rscript -e "library(rmarkdown);library(tidyverse);library(plotly)"
if [ $? -ne 0 ]
then
  echo "Some R packages are missing check rmarkdown, tidyverse and plotly are installed."
  exit 1
fi

# Get exe path for hicup
exePathForR=$(which R)
exePathForB2=$(which bowtie2)

# Define the output of HiCUP Digester
pathForDigest="${dirPathWithResults}/${genome}_digester_${restName}.txt.gz"


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

# Get the sample name and fastq file from the table
sample=$(cat ${filePathForTable} | awk -v i=${SLURM_ARRAY_TASK_ID} 'NR==i{print $1}')
relFilePathFastqR1=$(cat ${filePathForTable} | awk -v i=${SLURM_ARRAY_TASK_ID} 'NR==i{print $2}')
relFilePathFastqR2=$(cat ${filePathForTable} | awk -v i=${SLURM_ARRAY_TASK_ID} 'NR==i{print $3}')

# Each sample is processed into an independent directory:
pathResults=${dirPathWithResults}/${sample}/

# The directory is created (if not existing)
mkdir -p ${pathResults}

# The name of the sample is written in stdout
echo ${sample}

# The analysis part takes part within the pathResults
cd ${pathResults}

# Check if an output bam exists
inputBAM=$(find . -name "*.hicup.bam")

# We don't perform adapter removal step
# If there is no restriction site before adapter
# The pair will not be a valid pair

# Only run hicup if the bam does not exists
if [ -z $inputBAM ]; then
  if [ ! -e ${dirPathForFastq}/${relFilePathFastqR1} ]; then
    # If the fastq does not exists we assume it was an SRA ID
    mkdir -p ${dirPathForFastq}
    cd ${dirPathForFastq}
    # Write version to stdout:
    fasterq-dump --version
    if [ $? -ne 0 ]
    then
      echo "fasterq-dump is not installed and fastqFile not found so assumed it was a SRA ID.
Please install it for example in the conda environment (sra-tools>=2.11)."
      exit 1
    fi
    fasterq-dump -o ${sample}.fastq ${relFilePathFastqR1}
    if [ ! -s ${sample}_1.fastq ]; then
        echo "FASTQ R1 IS EMPTY"
        exit 1
    fi
    gzip ${sample}_1.fastq
    gzip ${sample}_2.fastq
    cd $pathResults
    relFilePathFastqR1=${sample}_1.fastq.gz
    relFilePathFastqR2=${sample}_2.fastq.gz
  fi
  if [ ! -s ${dirPathForFastq}/${relFilePathFastqR1} ]; then
    echo "FASTQ R1 IS EMPTY"
    exit 1
  fi
  # Run hicup
  hicup ${optionForHiCUP} --bowtie2 ${exePathForB2} \
    --digest ${pathForDigest} --format Sanger --index ${basenamePathForB2Index} \
    --keep --threads ${nbOfThreads} --zip --r ${exePathForR} \
    ${dirPathForFastq}/${relFilePathFastqR1} \
    ${dirPathForFastq}/${relFilePathFastqR2}
  # Update the inputBAM variable
  inputBAM=$(find . -name "*.hicup.bam")
else
  echo "HiCUP bam already exists. Not regenerating it."
fi

# Copy all final files to specific directories
mkdir -p ${dirPathWithResults}/allFinalFiles/reports
cp *.html ${dirPathWithResults}/allFinalFiles/reports/
