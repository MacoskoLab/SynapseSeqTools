import os
import argparse


# parse args
args_parser = argparse.ArgumentParser()
args_parser.add_argument("-i", "--umi_collapse_file", help="File with commands to run")
args = args_parser.parse_args()

umi_collapse_file = args.umi_collapse_file

# read in commands
cmds = []
with open(umi_collapse_file, "r") as f:
    for line in f:
        cmds.append(line.strip())

for umi_collapse_cmd in cmds:
    os.system(umi_collapse_cmd)