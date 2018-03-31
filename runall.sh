##################################
script_folder='/Users/gzx103/Documents/zhang_lab/projects/scripts/snapshot/bin/'
input_folder='/Users/gzx103/Documents/zhang_lab/projects/scripts/snapshot/test_data/input_data/'
output_folder='/Users/gzx103/Documents/zhang_lab/projects/scripts/snapshot/test_data/output_result/'


### run snapshot (CORE!!!)
echo 'run snapshot :o'
time python $script_folder'snapshot.py' -p peak_list.txt -n atac_4cell -t 5 -s signal_list.txt -l F -z F -x 0.01 -f function_list.txt -m mostfreq -c function_color_list.txt -e cd_tree.txt -i $input_folder -o $output_folder -b $script_folder
echo 'complete :)'
