#This script will create a folder named by the sample to create temporary files and will put the fastq.gz in fromGEO folder with the name of the sample_R1.fastq.gz and _R2.fastq.gz
pathForSRATable=$1
index=$2

SRAnumbers=(`cat $pathForSRATable | awk -v i=$index 'NR==i{print $3}'|tr "," " "`)
SlotNB=(`cat $pathForSRATable | awk -v i=$index 'NR==i{print $4}'|tr "," " "`)
sample=`cat $pathForSRATable | awk -v i=$index 'NR==i{print $2}'`
totNB=`echo ${SlotNB[@]} | awk '{for(i=1;i<=NF;i++){S+=$i}}END{print S}'`
if [ -e fromGEO/${sample}_R1.fastq.gz ]; then
  nbReads=`zcat fromGEO/${sample}_R1.fastq.gz | wc -l | awk '{print $1/4}'`
  if [ $nbReads = $totNB ]; then
    echo "fastq R1 already correct"
    if [ -e fromGEO/${sample}_R2.fastq.gz ]; then
      nbReads=`zcat fromGEO/${sample}_R2.fastq.gz | wc -l | awk '{print $1/4}'`
      if [ $nbReads = $totNB ]; then
        echo "fastq R2 already correct"
        exit 0
      else
        echo "fastq R2 incorrect"
        rm fromGEO/${sample}_R1.fastq.gz
        rm fromGEO/${sample}_R2.fastq.gz
      fi
    fi
  else
    echo "fastq R1 incorrect"
    rm fromGEO/${sample}_R1.fastq.gz
    if [ -e fromGEO/${sample}_R2.fastq.gz ]; then
      rm fromGEO/${sample}_R2.fastq.gz
    fi
  fi
fi
mkdir -p tmp_$sample
if [ -e tmp_${sample}/failed ]; then
 rm tmp_${sample}/failed
fi
nbSRA=${#SRAnumbers[@]}
for ((i=0; i<$nbSRA;i++))
do
  echo "sraName=${SRAnumbers[$i]}">tmp_${sample}/bash_${i}.sh
  echo "nbSpots=${SlotNB[$i]}">>tmp_${sample}/bash_${i}.sh
  echo "if [ -e tmp_${sample}/\${sraName}_1.fastq.gz ]; then">>tmp_${sample}/bash_${i}.sh
  echo "  # We first check if the existing R1 fastq is full">>tmp_${sample}/bash_${i}.sh

  echo "  nbSpotsWritten=\`zcat tmp_${sample}/\${sraName}_1.fastq.gz | wc -l | awk '{print \$1/4}'\`">>tmp_${sample}/bash_${i}.sh
  echo "  if [ ! \$nbSpots = \$nbSpotsWritten ]; then">>tmp_${sample}/bash_${i}.sh
  echo "    echo \"The file R1 present for \$sraName is not correct\"">>tmp_${sample}/bash_${i}.sh
  echo "    rm tmp_${sample}/\${sraName}*">>tmp_${sample}/bash_${i}.sh
  echo "  else">>tmp_${sample}/bash_${i}.sh
  echo "    echo \"The file R1 present for \$sraName is correct\"">>tmp_${sample}/bash_${i}.sh
  echo "    if [ -e tmp_${sample}/\${sraName}_2.fastq.gz ]; then">>tmp_${sample}/bash_${i}.sh
  echo "      nbSpotsWritten=\`zcat tmp_${sample}/\${sraName}_2.fastq.gz | wc -l | awk '{print \$1/4}'\`">>tmp_${sample}/bash_${i}.sh
  echo "      if [ ! \$nbSpots = \$nbSpotsWritten ]; then">>tmp_${sample}/bash_${i}.sh
  echo "        echo \"The file R2 present for \$sraName is not correct\"">>tmp_${sample}/bash_${i}.sh
  echo "        rm tmp_${sample}/\${sraName}*">>tmp_${sample}/bash_${i}.sh
  echo "      else">>tmp_${sample}/bash_${i}.sh
  echo "        echo \"The file R2 present for \$sraName is correct\"">>tmp_${sample}/bash_${i}.sh
  echo "        exit 0">>tmp_${sample}/bash_${i}.sh
  echo "      fi">>tmp_${sample}/bash_${i}.sh
  echo "    fi">>tmp_${sample}/bash_${i}.sh
  echo "  fi">>tmp_${sample}/bash_${i}.sh
  echo "fi">>tmp_${sample}/bash_${i}.sh
  echo "  # The quickest way to get fastq from sra is to use fasterq-dump.">>tmp_${sample}/bash_${i}.sh
  echo "echo \"Will use fasterq-dump for \$sraName\"">>tmp_${sample}/bash_${i}.sh
  echo "fasterq-dump \$sraName -O tmp_${sample}">>tmp_${sample}/bash_${i}.sh
  echo "if [ -e tmp_${sample}/\${sraName}_1.fastq ] && [ -s tmp_${sample}/\${sraName}_1.fastq ]; then">>tmp_${sample}/bash_${i}.sh
  echo "  gzip tmp_${sample}/\${sraName}_1.fastq &">>tmp_${sample}/bash_${i}.sh
  echo "  gzip tmp_${sample}/\${sraName}_2.fastq &">>tmp_${sample}/bash_${i}.sh
  echo "  wait">>tmp_${sample}/bash_${i}.sh
  echo "else">>tmp_${sample}/bash_${i}.sh
  echo "  echo \"fasterq failed will use fastq.\"">>tmp_${sample}/bash_${i}.sh
  echo "  # Else, we need to extract on the flow">>tmp_${sample}/bash_${i}.sh
  echo "  # Need to change home because fastq-dump will download in home/ncbi/public... and then there is no more place in home.">>tmp_${sample}/bash_${i}.sh
  echo "  export HOME=`pwd`">>tmp_${sample}/bash_${i}.sh
  echo "  fastq-dump --split-files --log-level fatal --accession \${sraName} --defline-seq '@\$sn[_\$rn]/\$ri' --defline-qual '+' --gzip -O tmp_${sample}">>tmp_${sample}/bash_${i}.sh
  echo "fi">>tmp_${sample}/bash_${i}.sh
  echo "nbSpotsWritten=\`zcat tmp_${sample}/\${sraName}_1.fastq.gz | wc -l | awk '{print \$1/4}'\`">>tmp_${sample}/bash_${i}.sh
  echo "if [ ! \$nbSpots = \$nbSpotsWritten ]; then">>tmp_${sample}/bash_${i}.sh
  echo "  echo \"The number of spots expected:${nbSpots} is different from what is in the fastq:${nbSpotsWritten}\"">>tmp_${sample}/bash_${i}.sh
  echo "  touch tmp_${sample}/failed">>tmp_${sample}/bash_${i}.sh
  echo "fi">>tmp_${sample}/bash_${i}.sh
  bash tmp_${sample}/bash_${i}.sh &
done
wait
if [ -e tmp_${sample}/failed ]; then
  exit 1
fi
cat tmp_${sample}/*_1.fastq.gz > fromGEO/${sample}_R1.fastq.gz
cat tmp_${sample}/*_2.fastq.gz > fromGEO/${sample}_R2.fastq.gz
rm -r tmp_${sample}
