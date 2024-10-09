import os
import argparse
import pickle
import gzip
import tqdm
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import scipy.io
import multiprocessing as mp
import itertools
from collections import Counter, defaultdict
import sys
sys.path.append("/home/jsilverm/06_synapseseq_repo/synapse_seq_pipeline_code/09_dedup_aav_lib")
import synapse_seq_functions as ssf

def write_mock_fastq(group_name, umi_vts, size_umi, size_vt, mock_fastq_folder):
    current_fastq_path = os.path.join(mock_fastq_folder, f"{group_name}.fastq")
    # create mock fastq file
    read_size = size_umi + size_vt
    # Write mock quality string for formatting, but is not taken into account for algo
    quality_string = "j" * read_size
    with open(current_fastq_path, "w") as fh:
        for umi_vt in umi_vts:
            umi = umi_vt["umi"]
            vt = umi_vt["vt"]
            count = umi_vt["count"]
            umi_vt = umi + vt
            for i in range(count):
                fh.write(f"@{group_name}\n")
                fh.write(f"{umi_vt}\n")
                fh.write(f"+\n")
                fh.write(f"{quality_string}\n")
    return current_fastq_path

def run_umi_vt_collapse(group_name, mock_fastq_folder, output_dir, log_base = "/home/jsilverm/logs"):
    mock_fastq_in = os.path.join(mock_fastq_folder, f"{group_name}.fastq")
    assert os.path.exists(mock_fastq_in)
    mock_fastq_out = os.path.join(output_dir, f"{group_name}_out.fastq")
    log_file = os.path.join(log_base, f"{group_name}_dedup_graph.log")
    cmd = f"umicollapse fastq -k 1 -i {mock_fastq_in} -o {mock_fastq_out} --tag > {log_file} 2>&1"
    os.system(cmd)


parser = argparse.ArgumentParser()
parser.add_argument("-b", "--base_dir", help="", required=True)
parser.add_argument("-r", "--region_name", help="", required=True)
parser.add_argument("-o", "--out_dir_general", help="", required=True)
parser.add_argument("-u", "--run_id", help="", required=True)
parser.add_argument("-t", "--is_test", action="store_true", required=False, default=False)
args = parser.parse_args()

base_dir = args.base_dir
region_name = args.region_name
out_dir_general = args.out_dir_general
run_id = args.run_id
is_test = args.is_test


print(f"base_dir: {base_dir}")
print(f"region_name: {region_name}")
print(f"out_dir_general: {out_dir_general}")
print(f"run_id: {run_id}")
print(f"is_test: {is_test}")



region_in_path = os.path.join(base_dir, region_name)
print(f"Checking that {region_in_path} exists")
assert os.path.exists(region_in_path), f"File {region_in_path} does not exist"

# Create output dir if it doesn't exist already
out_dir = os.path.join(out_dir_general, region_name, run_id)
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

print(f"Working out_dir: {out_dir}")

intermediate_dir = os.path.join(out_dir,"intermediate")
if not os.path.exists(intermediate_dir):
    os.mkdir(intermediate_dir)

# Create intermediate dir if it doesn't exist already. Used to store computed files
intermediate_dir_path = os.path.join(out_dir, "intermediate")
if not os.path.exists(intermediate_dir_path):
    os.makedirs(intermediate_dir_path)

sys.stderr.flush()
sys.stdout.flush()

# Hard coded local file names
reads_per_umi_fname = "raw_read_umi_matrix.mtx.gz"
vt_listfname = "raw_tags.txt.gz"
umi_listfname = "raw_umis.txt.gz"

# Read in umi-vt matrix to be umixvt in csr form
reads_per_umi_fname_full = os.path.join(region_in_path, reads_per_umi_fname)
print(f"Reading {reads_per_umi_fname_full}")
reads_per_umi_mat = scipy.io.mmread(reads_per_umi_fname_full)
reads_per_umi_mat = reads_per_umi_mat.tocsr()
umi_vt_mat = reads_per_umi_mat.T.tocsr()
print(reads_per_umi_mat.shape)

# Read in list of VTs and UMIs
vt_listfname_full = os.path.join(region_in_path, vt_listfname)
print(f"Reading {vt_listfname_full}")
vts = []
with gzip.open(vt_listfname_full, "rt") as f:
    for line in f:
        vts.append(line.strip())

umi_listfname_full = os.path.join(region_in_path, umi_listfname)
print(f"Reading {umi_listfname_full}")
umis = []
with gzip.open(umi_listfname_full, "rt") as f:
    for line in f:
        umis.append(line.strip())

bulk_obj = {"umi_vt_mat": umi_vt_mat, "umis": umis, "vts": vts}


sys.stderr.flush()
sys.stdout.flush()


# If test subset matrix to first 1000 umi and vts
if is_test:
    print("IS TEST. SUbsetting UMI and VTs")
    bulk_obj["umi_vt_mat"] = bulk_obj["umi_vt_mat"][:1000, :1000]
    bulk_obj["umis"] = bulk_obj["umis"][:1000]
    bulk_obj["vts"] = bulk_obj["vts"][:1000]

# Create mock fastq and output dedup folders
mock_fastq_folder = os.path.join(intermediate_dir, "mock_fastq")
if not os.path.exists(mock_fastq_folder):
    os.makedirs(mock_fastq_folder)

umi_dedup_outdir = os.path.join(intermediate_dir, "umi_dedup_fastqs")
if not os.path.exists(umi_dedup_outdir):
    os.makedirs(umi_dedup_outdir)



# Running hamming correction on combo of umi and corrected vt
# Create list of umi-vt pairs and # reads to write to a mock fastq
bulk_fastq_list = []
n_rows_umi_mat = vt_dedup_obj["umi_vt_mat"].shape[0]
for i in range(n_rows_umi_mat):
    umi = vt_dedup_obj["umis"][i]
    nonzero_elements = vt_dedup_obj["umi_vt_mat"][i, :].nonzero()[1]
    for nonzero_element in nonzero_elements:
        vt = vt_dedup_obj["vts"][nonzero_element]
        read_count = vt_dedup_obj["umi_vt_mat"][i, nonzero_element]
        umi_vt = {"umi": umi, "vt": vt, "count": read_count}
        bulk_fastq_list.append(umi_vt)


# Create mock fastqs for umicollapse in parallel
print("Writing mock fastq")
write_mock_fastq(region_name, bulk_fastq_list, umi_len,vt_len, mock_fastq_folder)

# Run umicollapse
print("Running umicollapse")
run_umi_vt_collapse(region_name, mock_fastq_folder, umi_dedup_outdir)

sys.stderr.flush()
sys.stdout.flush()

print("Parsing umicollapse results")
# Parse umicollapse results
# bulk_umi_dedup_outs will be a dictionary mapping observed umi-vt to corrected umi-vt
bulk_umi_dedup_outs = ssf.parse_umi_collapse_results(region_name, umi_dedup_outdir)


# Create a new u_mat with the corrected UMI sequences
# Use size of umi and vt to split the umi-vt string
u_mat_dedup_vt_umi_data_dict = defaultdict(Counter)
for i in tqdm.tqdm(range(vt_dedup_obj["umi_vt_mat"].shape[0])):
    current_umi = vt_dedup_obj["umis"][i]

    nonzero_vals = vt_dedup_obj["umi_vt_mat"][i, :].nonzero()[1]
    for nonzero_val in nonzero_vals:
        vt = vt_dedup_obj["vts"][nonzero_val]
        count = vt_dedup_obj["umi_vt_mat"][i, nonzero_val]

        umi_vt = current_umi + vt

        umi_vt_correction = bulk_umi_dedup_outs.get(umi_vt, None)
        if umi_vt_correction is None:
            print(f"None for {umi_vt}")
        assert umi_vt_correction is not None
        umi_corrected = umi_vt_correction[:umi_len]
        vt_corrected = umi_vt_correction[umi_len:]

        u_mat_dedup_vt_umi_data_dict[umi_corrected][vt_corrected] += count

sys.stderr.flush()
sys.stdout.flush()

# Create a new u_mat with the corrected UMI sequences
umis_from_dict = sorted(list(u_mat_dedup_vt_umi_data_dict.keys()))
vts_from_dict = sorted(set(i for v in u_mat_dedup_vt_umi_data_dict.values() for i in v.keys()))
u_mat_umi_dedup = ssf.create_matrix(u_mat_dedup_vt_umi_data_dict, umis_from_dict, vts_from_dict)
u_mat_vt_and_umi_dedup_obj = {"u": u_mat_umi_dedup, "umis": umis_from_dict, "vts": vts_from_dict}

# For each umi only keep reads for max vt
u_mat_dedup_max = ssf.filter_to_max_vt_per_cbumi(u_mat_vt_and_umi_dedup_obj["u"], u_mat_vt_and_umi_dedup_obj["umis"], u_mat_vt_and_umi_dedup_obj["vts"])
u_mat_dedup_max_obj = {"u": u_mat_dedup_max, "umis": u_mat_vt_and_umi_dedup_obj["umis"], "vts": u_mat_vt_and_umi_dedup_obj["vts"]}

# Save pre read filt object
u_mat_out = os.path.join(intermediate_dir, "u_mat_dedup_max.pkl")
print(f"Writing u_mat to {u_mat_out}")
with open(u_mat_out, "wb") as fh:
    pickle.dump(u_mat_dedup_max_obj, fh)


# Create plots and save tags for different read and umi filters
pdf_out_name = os.path.join(out_dir, f"{region_name}_dedup_summary.pdf")
pp = PdfPages(pdf_out_name)

reads_umi_fig, ax = plt.subplots()
reads_per_umi = ssf.create_reads_umi_histogram(u_mat_dedup_max_obj["u"], title_base = f"Bulk {region_name}", truncate_value=5000, bins=80, ax=ax)
pp.savefig(reads_umi_fig)
plt.close(reads_umi_fig)


# remove the mock fastq created and gzip the deduped fastqs
print("removing mock fastq folder")
os.system(f"rm -r {mock_fastq_folder}")
print("gzipping deduped fastqs")
os.system(f"find {umi_dedup_outdir} -type f -exec gzip {{}} \;")

print("Finished")



