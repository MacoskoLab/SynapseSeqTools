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


# get dirname and create dir
os.makedirs(os.path.dirname(out_path), exist_ok=True)

# get current file path
current_dir = os.path.dirname(os.path.abspath(__file__))

# get one dir above
base_ss_dir = os.path.dirname(current_dir)
config_yaml_path = os.path.join(base_ss_dir, "config.yaml")
config = yaml.safe_load(open(config_yaml_path, "r"))

bulk_parse_fastq_path = os.path.join(base_ss_dir, "02_bulk", "extract_barcodes_bulk.py")
assert os.path.exists(bulk_parse_fastq_path), f"bulk_parse_fastq_path {bulk_parse_fastq_path} does not exist"

prepend_path = config["prepend_cmds_path"]

# read in samplesheet
sample_df = pd.read_csv(args.samplesheet, sep=',')

cmds_list = []

for inx, row in sample_df.iterrows():
    fastq_dir = row['fastq_dir']
    # read full path for r1 and r2
    # list files in fastq_dir and filter for r1 and r2
    r1 = [f for f in os.listdir(fastq_dir) if "_R1_" in f][0]
    r2 = [f for f in os.listdir(fastq_dir) if "_R2_" in f][0]
    r1_path = os.path.join(fastq_dir, r1)
    r2_path = os.path.join(fastq_dir, r2)

    assert os.path.exists(r1_path), f"r1 file {r1_path} does not exist"
    assert os.path.exists(r2_path), f"r2 file {r2_path} does not exist"

    out_dir_extract = row['out_dir_extract']
    sample_id = row['sample_id']
    out_dir = os.path.join(out_dir_extract, sample_id)
    os.makedirs(out_dir, exist_ok=True)

    # create command
    out_log_file = os.path.join(out_dir, f"{sample_id}.log")
    cmd = f"python {bulk_parse_fastq_path} {r1_path} {r2_path} --output-dir {out_dir} > {out_log_file} 2>&1"
    if len(prepend_path) > 0:
        cmd = f"{prepend_path} '{cmd}'"
    cmd_tup = (sample_id, cmd)
    cmds_list.append(cmd_tup)


if use_slurm:
    for sample_id, cmd in cmds_list:
        # create sbatch file
        sbatch_path = os.path.join(out_path, f"{sample_id}.sh")

        out_path = os.path.join(out_path, f"{sample_id}.out")
        err_path = os.path.join(out_path, f"{sample_id}.err")
        MEM = "100G"
        partition = "hpcx_macosko"
        with open(sbatch_path, "w") as f:
            f.write(f"#!/bin/bash\n")
            f.write(f"#SBATCH --job-name={sample_id}\n")
            f.write(f"#SBATCH --output={out_path}\n")
            f.write(f"#SBATCH --error={err_path}\n")
            f.write(f"#SBATCH --time=24:00:00\n")
            f.write(f"#SBATCH --cpus-per-task=2\n")
            f.write(f"#SBATCH --mem={MEM}\n")
            f.write(f"#SBATCH --nodes=1\n")
            f.write(f"#SBATCH --ntasks=1\n")
            if partition is not None:
                f.write(f"#SBATCH --partition={partition}\n")
            f.write(f"{cmd}\n")
        # run sbatch file
        # subprocess.run(f"sbatch {sbatch_path}", shell=True)
else:
    # write commands to file
    with open(out_path, "w") as f:
        for cmd in cmds_list:
            f.write(f"{cmd[1]}\n")


