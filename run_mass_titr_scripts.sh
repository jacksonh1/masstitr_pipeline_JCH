# ==============================================================================
# // parameters
# ==============================================================================


#== multiplexed forward and reverse reads ==#
forward_read_file="./data/sample_data_forward.fastq"
reverse_read_file="./data/sample_data_reverse.fastq"

#== directory to output unsorted barcode files after de-multiplexing the forward and reverse fastq files ==#
barcode_directory="./data/"

#== directory to put folders with sorted barcode files. The script creates separate folders for each `experiment` as defined in `sample_key_file` and moves the barcode files in `barcode_directory` into the correct `experiment` folder ==#
experiment_directory="./"

sample_key_file="FINAL_sample_key_cell_counts_allreps.csv"

#== bbmap `reformat` executable ==#
reformat_command="reformat.sh"

#== output folder to save binding curves ==#
# binding_curve_folder="./binding_curves_weighted_3_param_read_count_longform"
# binding_curve_folder2="./binding_curves_weighted_3_param_read_count_AAclustering_longform"
binding_curve_folder="./binding_curves_weighted_3_param_cell_count_AAclustering_longform"
# ==============================================================================
# // Run scripts
# ==============================================================================

~/anaconda3/envs/py2/bin/python ./src/mt01_barcode_sort.py $forward_read_file $reverse_read_file $barcode_directory
~/anaconda3/envs/py2/bin/python ./src/mt02_sort_barcodes.py $barcode_directory $experiment_directory $sample_key_file
~/anaconda3/envs/py2/bin/python ./src/mt03_process_barcode_files.py $experiment_directory $sample_key_file $reformat_command


~/anaconda3/envs/main/bin/python ./src/mt04_process_titrations_new.py $experiment_directory $sample_key_file $binding_curve_folder
~/anaconda3/envs/main/bin/python ./src/mt05_fit_kd_longform.py "$binding_curve_folder/exp1-binding_curves.csv"
~/anaconda3/envs/main/bin/python ./src/mt05_fit_kd_longform.py "$binding_curve_folder/exp2-binding_curves.csv"
~/anaconda3/envs/main/bin/python ./src/mt05_fit_kd_longform.py "$binding_curve_folder/exp3-binding_curves.csv"
~/anaconda3/envs/main/bin/python ./src/mt05_fit_kd_longform.py "$binding_curve_folder/exp4-binding_curves.csv"
# ~/anaconda3/envs/main/bin/python ./src/correlation_heatmaps_noR2.py "$binding_curve_folder"

