# snapshot


### Dependence:
#### Python/2.7
###### numpy; subprocess; collections; getopt; sys; os
#### R
###### ggplot2; pheatmap; igraph; networkD3

### input data
###### The cell type peak binary label file list: 1st column is the foldername and the filename in input folder; 2nd column is the cell type label in output figures
```
peak_list.txt
head peak_list.txt 
atac_pk/LSK.pk.bed	LSK
atac_pk/CMP.pk.bed	CMP
atac_pk/MEP.pk.bed	MEP
atac_pk/GMP.pk.bed	GMP
```

###### The cell type peak signal file list: 1st column is the foldername and the filename in input folder; 2nd column is the cell type label in output figures
```
signal_list.txt
head signal_list.txt 
atac_sig/LSK.atac.sig.bed	LSK
atac_sig/CMP.atac.sig.bed	CMP
atac_sig/MEP.atac.sig.bed	MEP
atac_sig/GMP.atac.sig.bed	GMP
```

###### The cell type functional state file list: 1st column is the foldername and the filename in input folder; 2nd column is the cell type label in output figures
```
function_list.txt
head function_list.txt 
ideas_label/LSK.ideas.bed	LSK
ideas_label/CMP.ideas.bed	CMP
ideas_label/MEP.ideas.bed	MEP
ideas_label/GMP.ideas.bed	GMP
```

###### The cell type differentiation tree: Each row represent one edge in the ell type differentiation tree. The 1st cell type is the progenitor cell type and the 2nd cell type is the differentiated cell type
```
cd_tree.txt 
head cd_tree.txt
LSK,CMP
CMP,MEP
CMP,GMP
```

###### The functional state color list
```
head function_color_list.txt
36	35	194,7,153	250,151,3
35	34	250,151,3	136,53,241
34	33	136,53,241	197,151,0
33	32	197,151,0	138,177,89
32	31	138,177,89	191,0,84
31	30	191,0,84	176,0,93
30	29	176,0,93	252,48,50
29	28	252,48,50	0,0,172
28	27	0,0,172	219,8,0
27	26	219,8,0	241,198,171
```

### snapshot
###### change the folder names in 'runall.sh'
```
script_folder='/snapshot_script_folder/snapshot/bin/'
input_folder='/snapshot_script_folder/snapshot/test_data/input_data/'
output_folder='/snapshot_script_folder/snapshot/test_data/output_result/'
```

###### run snapshot
```
bash runall.sh
```

