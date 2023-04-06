# Get the SraRunTable.txt for GSE96107
path=/home/ldelisle/softwares/Hi-C_reanalysis_Bonev_2017/scripts_mm39/
GSM_to_keep="18 19 20 21 35 36 37 38 39 40 41 42"
pathForTable=${path}/HiC_Bonev.txt
pathForTable2=${path}/HiC_Bonev_GSM.txt
for gsm in $GSM_to_keep; do
    awk -F "," -v i=$gsm 'BEGIN{j=1}$29=="GSM25338"i{print $29"_"j"\t"$1; j+=1}' ${path}/SraRunTable.txt > temp.txt
    echo -e "GSM25338${gsm}\t$(cut -f 1 temp.txt | tr "\n" ",")" >> ${pathForTable2}
    cat temp.txt >> ${pathForTable}
done
rm temp.txt
sbatch /home/ldelisle/softwares/Hi-C_reanalysis_Bonev_2017/scripts_mm39/1_HiC_hicup_on_each_sra.sh
sbatch /home/ldelisle/softwares/Hi-C_reanalysis_Bonev_2017/scripts_mm39/2_HiC_dedup_per_replicate_custom.sh
sbatch /home/ldelisle/softwares/Hi-C_reanalysis_Bonev_2017/scripts_mm39/3_HiC_merge_rep.sh
sbatch /home/ldelisle/softwares/Hi-C_reanalysis_Bonev_2017/scripts_mm39/4_HiC_generate_matrices.sh
# Submitted batch job 2786664
# I forgot the multiplecorrection param:
sbatch --dependency=afterany:2786819 /home/ldelisle/softwares/Hi-C_reanalysis_Bonev_2017/scripts_mm39/HiC_call_TADs.sh

sbatch /home/ldelisle/softwares/Hi-C_reanalysis_Bonev_2017/scripts_mm39/HiC_merge_rep_TADs.sh

sbatch /home/ldelisle/softwares/Hi-C_reanalysis_Bonev_2017/scripts_mm39/HiC_call_TADs_500.sh

sbatch /home/ldelisle/softwares/Hi-C_reanalysis_Bonev_2017/scripts_mm39/HiC_merge_rep_TADs_500.sh

sbatch /home/ldelisle/softwares/Hi-C_reanalysis_Bonev_2017/scripts_mm39/HiC_call_TADs_20kb.sh
sbatch /home/ldelisle/softwares/Hi-C_reanalysis_Bonev_2017/scripts_mm39/HiC_merge_rep_TADs_20kb.sh 

sbatch /home/ldelisle/softwares/Hi-C_reanalysis_Bonev_2017/scripts_mm39/HiC_call_TADs_20kb_500.sh 
sbatch --dependency=afterany:2806622 /home/ldelisle/softwares/Hi-C_reanalysis_Bonev_2017/scripts_mm39/HiC_merge_rep_TADs_20kb_500.sh

sbatch /home/ldelisle/softwares/Hi-C_reanalysis_Bonev_2017/scripts_mm39/HiC_call_TADs_20kb_800.sh 
sbatch --dependency=afterany:2806849 /home/ldelisle/softwares/Hi-C_reanalysis_Bonev_2017/scripts_mm39/HiC_merge_rep_TADs_20kb_800.sh