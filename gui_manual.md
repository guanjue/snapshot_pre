# Snapshot-GUI User Manual

## Prerequisites
- pyqt5
Snapshot-GUI is based on PyQt5. You can install PyQt5 using command `pip install pyqt5`, or install it from [PyQt5 official website](https://www.riverbankcomputing.com/software/pyqt/download5)

## Get Start

To start Snapshot-GUI, run command ```python ./bin/gui/main.py```.
Figure 1 illustrate the Snapshot-GUI start window.

#### Figure 1
<img src="https://github.com/guanjue/snapshot/blob/master/test_data/example/gui/start.png" width="800"/>

## Use Snapshot-GUI

### 1. choose input directory.

Click the `Choose Directory` button for `Input Directory` to choose input directory.

The input directory should contains necessary input files and directories, as illustrated in Figure 2:

- `function_label` fold
- `atac_sig` fold
- `atac_pick` fold
- `signal_list.txt` file
- `peak_list.txt` file
- `function_list.txt` file
- `function_color_list.txt` file
- `cd_tree.txt` file

#### Figure 2
<img src="https://github.com/guanjue/snapshot/blob/master/test_data/example/gui/choose_input.png" width="800"/>

### 2. choose output directory.

Click the `Choose Directory` button for `Output Directory` to choose output directory. All analysis results will be saved in output directory.


### 3. preprocessing

Click the `Preprocess` button to check the input and output directories. Please make sure the input and output directories are chosen before preprocessing, otherwise, you will get warning as in Figure 3. 

#### Figure 3
<img src="https://github.com/guanjue/snapshot/blob/master/test_data/example/gui/preprocess.png" width="800"/>

After preprocessing, the `OK` button is able to use, and the `Change Display Color` and `Change Parameters` buttons are also able to use, as illustrated in Figure 4. 

#### Figure 4
<img src="https://github.com/guanjue/snapshot/blob/master/test_data/example/gui/ready_to_run.png" width="800"/>

### 4. change parameters and display colors

You can change parameters and display colors by click `Change Parameters` and `Change Display Color` button.

### 5. analyze your data

Click `OK` button to analyze your data, the Snapshot-GUI will automatically analyze your data in background, the progress bar is used to indicate the analysis progress, as shown in Figure 5. 

#### Figure 5
<img src="https://github.com/guanjue/snapshot/blob/master/test_data/example/gui/processing.png" width="800"/>

### 6. visualize the result

When the data analysis finishes, the visualization window will pop out automatically, as illustrated in Figure 6.

#### Figure 6
<img src="https://github.com/guanjue/snapshot/blob/master/test_data/example/gui/visual_result.png" width="800"/>

The default visualization is Index-set mean signal heatmap (as in Figure 6). You can also choose to visualize the Index-set most frequent functional annotation heatmap, as shown in Figure 7.

#### Figure 7
<img src="https://github.com/guanjue/snapshot/blob/master/test_data/example/gui/change_visual.png" width="800"/>

### 7. collect results

All analysis results are stored in output directory, as illustrated in Figure 8. 

#### Figure 8
<img src="https://github.com/guanjue/snapshot/blob/master/test_data/example/gui/result_dir.png" width="800"/>
