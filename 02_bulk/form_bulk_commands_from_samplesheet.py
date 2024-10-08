import argparse
import os


parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="Input samplesheet file", required=True)
parser.add_argument("-o", "--output", help="Output bulk commands file", required=True)
args = parser.parse_args()

input_file = args.input
output_file = args.output

with open(input_file, "r") as f:
    col_to_value_dict = {}
    column_names = f.readline().strip().split(",")

    for col in column_names:
        col_to_value_dict[col] = []
    
    for row in f:
        row_entries = row.strip().split(",")
        for col, value in zip(column_names, row_entries):
            col_to_value_dict[col].append(value)


exec_file="/home/jsilverm/06_synapseseq_repo/synapse_seq_pipeline_code/09_dedup_aav_lib/bulk/bulk_dedup_full.py"
# form commands
cmds = []
for i in range(len(col_to_value_dict['base_dir'])):
    cmd = f"python {exec_file} -b {col_to_value_dict['base_dir'][i]} -r {col_to_value_dict['region_name'][i]} -o {col_to_value_dict['out_dir_general'][i]} -u {col_to_value_dict['run_id'][i]} > /home/jsilverm/logs/{col_to_value_dict['region_name'][i]}_bulk.log 2>&1"
    cmds.append(cmd)

# write commands to file
with open(output_file, "w") as f:
    for cmd in cmds:
        f.write(cmd + "\n")