import argparse
import os
import pickle
import subprocess
import datetime

# read args from command line
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_path", help="Input Path")
parser.add_argument("-s", "--input_sample", help="Input Sample")
parser.add_argument("-o", "--output_folder", help="Output Folder, will be made within input path")
parser.add_argument("-p", "--polygon_list", help="Polygon list. Sep by ;")

args = parser.parse_args()

input_path = args.input_path
input_sample = args.input_sample
output_folder = args.output_folder
polygon_list = args.polygon_list

# Will run 2 scripts, one for dedup and one for regional matching

expected_out_pkl_file = os.path.join(input_path, input_sample, output_folder, "dedup_vt_and_umi_obj.pkl")
print(f"Searching for dedup pkl file: {expected_out_pkl_file}")
if not os.path.exists(expected_out_pkl_file):
    print("Could not find dedup pkl file, running dedup script")
    # Dedup
    dedup_script_path="/home/jsilverm/06_synapseseq_repo/synapse_seq_pipeline_code/09_dedup_aav_lib/slideseq/slide_seq_hamming_dedup.py"
    dedup_cmd=f"python {dedup_script_path} -i {input_path} -s {input_sample} -o {output_folder} > /home/jsilverm/logs/{input_sample}_dedup.log 2>&1"
    print(dedup_cmd)
    subprocess.run(dedup_cmd, shell=True)
    # get success status

    if not os.path.exists(expected_out_pkl_file):
        raise Exception("Dedup script failed. Exiting")
else:
    print(f"Found dedup pkl file: {expected_out_pkl_file}")
    print("Skipping dedup step")

# match and read filter
working_dir = os.path.join(input_path, input_sample)
match_script_path="/home/jsilverm/06_synapseseq_repo/synapse_seq_pipeline_code/09_dedup_aav_lib/slideseq/slide_seq_region_matching.py"
filter_cmd = f"python {match_script_path} -w {working_dir} -b {output_folder} -l {polygon_list} > /home/jsilverm/logs/{input_sample}_region_match.log 2>&1"
print(filter_cmd)
subprocess.run(filter_cmd, shell=True)


