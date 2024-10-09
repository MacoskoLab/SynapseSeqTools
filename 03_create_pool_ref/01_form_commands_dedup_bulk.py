import pandas as pd
import numpy as np
import os
import subprocess
import logging
import argparse
import yaml

# read in arguments from command line

parser = argparse.ArgumentParser(description='Create pool reference')
parser.add_argument('--samplesheet', type=str, help='input file')
parser.add_argument('--output', type=str, help='output file')
# set as store true if passed
parser.add_argument('--use_slurm', action='store_true', help='use slurm')
args = parser.parse_args()

# IF using slurm, will write sbatch files to out_path, need to store this in a directory
use_slurm = args.use_slurm
out_path = args.output
# if use_slurm ensure out_path is a directory
if use_slurm:
    # assert out_path is not a file
    assert not os.path.isfile(out_path), f"out_path {out_path} is a file, should be a directory"
    # create directory if it does not exist
    os.makedirs(out_path, exist_ok=True)


# get current file path
current_dir = os.path.dirname(os.path.abspath(__file__))
# get one dir above
base_ss_dir = os.path.dirname(current_dir)
config_yaml_path = os.path.join(base_ss_dir, "config.yaml")
config = yaml.safe_load(open(config_yaml_path, "r"))

prepend_path = config["prepend_cmds_path"]

# read in samplesheet
sample_df = pd.read_csv(args.samplesheet, sep=',')

cmds_list = []
exec_path = os.path.join(current_dir, "02_graph_dedup_bulk_pool.py")

for inx, row in sample_df.iterrows():
    out_dir_extract = row['out_dir_extract']
    out_dir_dedup = row['out_dir_dedup']
    sample_id = row['sample_id']
    run_id = row['run_id']

    out_log_file = os.path.join(out_dir_dedup, f"{sample_id}_dedup.log")

    cmd = f"{prepend_path}python {exec_path} -b {out_dir_extract} -r {sample_id} -o {out_dir_dedup} -u {run_id} > {out_log_file} 2>&1"
    cmds_list.append(cmd)

# write to file
with open(out_path, "w") as fh:
    for cmd in cmds_list:
        fh.write(f"{cmd}\n")




