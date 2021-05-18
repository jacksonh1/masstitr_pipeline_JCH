barcode_output_directory="./data/"
forward_read_file="./data/sample_data_forward.fastq"
reverse_read_file="./data/sample_data_reverse.fastq"
experiment_directory="./"
sample_key_file="./FINAL_sample_key_w_cell_counts.csv"
reformat_command="reformat.sh"


python ./src/01_barcode_sort.py $forward_read_file $reverse_read_file $barcode_output_directory
python ./src/02_sort_barcodes.py $barcode_output_directory $experiment_directory $sample_key_file
python ./src/03_process_barcode_files.py $experiment_directory $sample_key_file $reformat_command
~/anaconda3/envs/main/bin/python ./src/04_process_titrations.py $experiment_directory $sample_key_file
