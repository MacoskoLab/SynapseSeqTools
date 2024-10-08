import os
import argparse
import pickle
import gzip



parser = argparse.ArgumentParser()
parser.add_argument("-b", "--base_dir", help="Input samplesheet file", required=True)
parser.add_argument("-r", "--region_name", help="Output bulk commands file", required=True)
parser.add_argument("-o", "--out_dir_general", help="Output bulk commands file", required=True)
parser.add_argument("-u", "--run_id", help="Output bulk commands file", required=True)
args = parser.parse_args()

base_dir = args.base_dir
region_name = args.region_name
out_dir_general = args.out_dir_general
run_id = args.run_id

print(f"base_dir: {base_dir}")
print(f"region_name: {region_name}")
print(f"out_dir_general: {out_dir_general}")
print(f"run_id: {run_id}")


region_in_path = os.path.join(base_dir, region_name)
print(f"Checking that {region_in_path} exists")
assert os.path.exists(region_in_path), f"File {region_in_path} does not exist"

out_dir = os.path.join(out_dir_general, region_name, run_id)
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

tags_fname = "raw_tags.txt.gz"
raw_umis_fname = "raw_umis.txt.gz"
raw_umi_tag_fname = "umis_per_tags.txt.gz"

vt_hamming_k = 1

intermediate_dir_path = os.path.join(out_dir, "intermediate")
if not os.path.exists(intermediate_dir_path):
    os.makedirs(intermediate_dir_path)

tags_fname_full = os.path.join(region_in_path, tags_fname)
raw_umis_fname_full = os.path.join(region_in_path, raw_umis_fname)
raw_umi_mat_fname_full = os.path.join(region_in_path, raw_umi_tag_fname)

assert os.path.exists(tags_fname_full), f"File {tags_fname_full} does not exist"
assert os.path.exists(raw_umis_fname_full), f"File {raw_umis_fname_full} does not exist"
assert os.path.exists(raw_umi_mat_fname_full), f"File {raw_umi_mat_fname_full} does not exist"

print(f"Reading {tags_fname_full}")
vts = []
with gzip.open(tags_fname_full, "rt") as f:
    for line in f:
        vts.append(line.strip())


vt_whitelist_path = "/home/jsilverm/06_synapseseq_repo/synapse_seq_pipeline_code/09_dedup_aav_lib/00_hamming_dist_objs/aav_whitelist_hamming2.txt"
with open(vt_whitelist_path, "r") as f:
    vt_whitelist = set([x.strip() for x in f.readlines()])

vts_in_whitelist = set([x for x in vts if x in vt_whitelist])

out_path = os.path.join(out_dir, "whitelist_vts.pkl")
with open(out_path, "wb") as f:
    pickle.dump(vts_in_whitelist, f)

