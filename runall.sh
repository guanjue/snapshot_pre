##################################
script_folder='/Users/universe/Documents/2018_BG/snapshot/bin/'

### run snapshot (CORE!!!)
echo 'run snapshot :o'
time python $script_folder'snapshot.py' -p peak_list.txt -t 200 -s signal_list.txt -l T -f ideas_list.txt -t str -c ideas_color_list.txt -i test_data -o test_data_output
echo 'complete :)'