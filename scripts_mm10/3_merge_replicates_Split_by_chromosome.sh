#!/bin/bash

#SBATCH -o slurm-%x-%A_%2a.out
#SBATCH -e slurm-%x-%A_%2a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lucille.delisle@epfl.ch
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --mem 10G
#SBATCH --cpus-per-task 28
#SBATCH --time 48:00:00
#SBATCH --array=1,33,89
#SBATCH --job-name SplitcsortvalidPairs
#SBATCH --workdir /scratch/ldelisle/HiCUP

pathForTableWithSRA="table.txt"
genome="mm10"
pathForSizes="/home/ldelisle/softwares/Hi-C_reanalysis_Bonev_2017/${genome}.size"

sampleOri=`cat $pathForTableWithSRA | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $2}'`
sampleBase=`basename $sampleOri _1_1`
name=${sampleBase}_Merge
mkdir -p $name

##NEEDED FOR BIG FILES##
mkdir -p temp
export TMPDIR=/scratch/ldelisle/HiCUP/temp

module purge
module load gcc/6.4.0 htslib/1.8

mkdir -p $name
cd $name
# First concatenate the valid pairs non sorted
if [ ! -e  ${name}.validPairs_notSorted.txt.gz ]; then
  for ((i=1;i<=4;i++)); do
    cat ../${sampleBase}_${i}/${sampleBase}_${i}.validPairs_nonSorted.txt.gz >> ${name}.validPairs_notSorted.txt.gz
  done
fi

# Split by chromosome for cis and put all trans together
nbChrGZ=`ls *_split_chr*validPairs_notSorted.txt.gz | wc -l`

if [ $nbChrGZ -eq 0 ]; then
  zcat ${name}.validPairs_notSorted.txt.gz | awk -v base="${name}_split" '$3==$7{print > base"_"$3".validPairs_notSorted.txt"}$3!=$7{print}' | gzip > ${name}_split_trans.validPairs_notSorted.txt.gz
  for f in ${name}_split_chr*.validPairs_notSorted.txt; do
    gzip -f $f &
  done
fi

wait

mkdir -p logs

# Use cooler csort to sort the pairs and index with tabix
# cooler version 0.7.11 was installed with pip

for f in ${name}_split_*.validPairs_notSorted.txt.gz; do
  tname=`basename $f .validPairs_notSorted.txt.gz`
  tpairs=${tname}.validPairs.csort.gz
  if [ ! -e ${tpairs}.tbi ]; then
    cooler csort -i tabix -c1 3 -c2 7 -p1 4 -p2 8 -o $tpairs $f $pathForSizes 2> logs/${tname}.log &
  fi
done

wait
# Put back the TMPDIR
export TMPDIR=/tmp
