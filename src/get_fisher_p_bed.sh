#cat /Volumes/MAC_Data/data/labs/zhang_lab/01projects/CD_viewer/test_data/ideas_label/20cell_nbnobp.LSK_BM.bed | awk -F '\t' -v OFS='\t' '{print $1,$2,$3, "200bins"}' > 20cell.bed
ls *.atacrep.fisher_p.txt > fisher_p_list.txt
for file in $(cat fisher_p_list.txt)
do
	paste 20cell.bed $file | awk -F '\t' -v OFS='\t' '{print $1,$2,$3,$4,$5,"+"}'> $file'.bed'
done
