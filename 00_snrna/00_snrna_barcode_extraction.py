import gzip
from pathlib import Path
import pickle
import numpy as np
import multiprocessing as mp
import os
import sys
import json
import argparse
sys.path.append("/home/jsilverm/06_synapseseq_repo/synapse_seq_pipeline_code/09_dedup_aav_lib/")

import synapse_seq_functions as ssf

def check_constant_seq_validity(genomic_seq, absolute_start_inx_wpre, 
                                absolute_end_inx_wpre, constant_sequence_wpre, 
                                absolute_start_inx_polya, absolute_end_inx_polya, 
                                constant_sequence_polya, wpre_hamming_limit=3, polya_hamming_limit=3):
    """
    For a given read (genomic_seq) check if the constant sequences are valid
    as defined by the indices for their start and end in the read as well as the hamming
    limits.
    Will be run in parallel.
    """
    current_wpre_pos_reads = genomic_seq[absolute_start_inx_wpre: absolute_end_inx_wpre]
    assert len(current_wpre_pos_reads) == len(constant_sequence_wpre)
    hamming_dist_wpre = ssf.do_basic_hamming(current_wpre_pos_reads, constant_sequence_wpre)
    valid_wpre = hamming_dist_wpre <= wpre_hamming_limit

    # check polya
    current_polya_pos_reads = genomic_seq[absolute_start_inx_polya: absolute_end_inx_polya]
    assert len(current_polya_pos_reads) == len(constant_sequence_polya)
    hamming_dist_polya = ssf.do_basic_hamming(current_polya_pos_reads, constant_sequence_polya)
    valid_polya = hamming_dist_polya <= polya_hamming_limit 

    is_valid_read = valid_wpre and valid_polya
    if is_valid_read:
        final_seq = genomic_seq
    else:
        final_seq = None

    res = {"is_valid_read": is_valid_read, "final_seq": final_seq, "hamming_dist_wpre": hamming_dist_wpre, "hamming_dist_polya": hamming_dist_polya}

    return res



# Work from merged fastq
# input_dir = "/mnt/disks/synapseseqfastq/V1/Mickey_A/"
# merged_fastq_name = "merged.fastq.gz"

#input dir also is where files are written
argparser = argparse.ArgumentParser()
argparser.add_argument("-i", "--input_dir", help="Directory with merged fastq files")
argparser.add_argument("-f", "--fastq_name", default="merged.fastq.gz" ,help="Name of merged fastq file")

args = argparser.parse_args()

input_dir = args.input_dir
merged_fastq_name = args.fastq_name
output_dir = input_dir

print("Input dir: ", input_dir)
print("Merged fastq name: ", merged_fastq_name)

sys.stderr.flush()
sys.stdout.flush()

# Define read structure variables
constant_sequence_wpre = "GATACCGAGCGCTGC"
constant_sequence_polya = "TCGAGAGATCTACGGG"

r1_length=50
r2_length=150

cellbc_start_r1=0
cellbc_length=16

umi_start_r1=16
umi_length=12

vt_r2_start=25
vt_length=32

wpre_const_r2_start=5
polya_const_r2_start=57

wpre_hamming_limit = 3
polya_hamming_limit = 3

n_cores=35

vt_f_name = "VTs.fastq.gz"
cb_f_name = "CBCs.fastq.gz"
umi_f_name = "UMIs.fastq.gz"

################################################
            # Load merged fastq #
################################################

merged_fastq_full_path = os.path.join(input_dir, merged_fastq_name)
if not os.path.exists(merged_fastq_full_path):
    print("Merged fastq not found at: ", merged_fastq_full_path)
    exit(1)

print(f"Loading merged fastq from {merged_fastq_full_path}")
with gzip.open(merged_fastq_full_path, 'rt') as f:
    # strip each line and remove \t
    genomic_seqs = [line.strip().replace("\t", "") for line_inx, line in enumerate(f)]


# define the start and end index of the constant sequences in R2
absolute_start_inx_wpre = r1_length + wpre_const_r2_start
absolute_end_inx_wpre = absolute_start_inx_wpre + len(constant_sequence_wpre)

absolute_start_inx_polya = r1_length + polya_const_r2_start
absolute_end_inx_polya = absolute_start_inx_polya + len(constant_sequence_polya)

print("absolute_start_inx_wpre: ", absolute_start_inx_wpre)
print("absolute_end_inx_wpre: ", absolute_end_inx_wpre)

print("absolute_start_inx_polya: ", absolute_start_inx_polya)
print("absolute_end_inx_polya: ", absolute_end_inx_polya)

sys.stderr.flush()
sys.stdout.flush()


#####################################################
# Read 2 Structure Filtering in Parallel #
#####################################################

print("Constructing args for parallel read processing")
read_processing_args = [(genomic_seq, absolute_start_inx_wpre, 
                         absolute_end_inx_wpre, constant_sequence_wpre, 
                         absolute_start_inx_polya, absolute_end_inx_polya, 
                         constant_sequence_polya) for genomic_seq in genomic_seqs]
print(len(read_processing_args))
pool = mp.Pool(n_cores)
print("Processing reads in parallel")
res = pool.starmap(check_constant_seq_validity, read_processing_args)

filtered_genomic_seqs = [r["final_seq"] for r in res if r["is_valid_read"]]
wpre_hamming_lis = [r["hamming_dist_wpre"] for r in res]
polya_hamming_lis = [r["hamming_dist_polya"] for r in res]


#####################################################
    # Extract Cell Barcodes, UMIs, and VTs #
#####################################################

# open the file and extract the cell barcodes, umis, and vt sequences
vts = []
cellbcs = []
umis = []

cb_end_inx = cellbc_start_r1 + cellbc_length
umi_end_inx = umi_start_r1 + umi_length

vt_absolute_start_inx = vt_r2_start + r1_length
vt_absolute_end_inx = vt_absolute_start_inx + vt_length

print(f"Reading CB from {cellbc_start_r1} to {cb_end_inx}")
print(f"Reading UMI from {umi_start_r1} to {umi_end_inx}")
print(f"Reading VT from {vt_absolute_start_inx} to {vt_absolute_end_inx}")

sys.stderr.flush()
sys.stdout.flush()

for genomic_seq in filtered_genomic_seqs:
    cellbc = genomic_seq[cellbc_start_r1: cb_end_inx]
    umi = genomic_seq[umi_start_r1: umi_end_inx]
    vt = genomic_seq[vt_absolute_start_inx: vt_absolute_end_inx]

    cellbcs.append(cellbc)
    umis.append(umi)
    vts.append(vt)


filter_stats = {
    "wpre_hamming_limit": wpre_hamming_limit,
    "polya_hamming_limit": polya_hamming_limit,
    "n_reads_start": len(genomic_seqs),
    "n_reads_end": len(filtered_genomic_seqs)
}

#####################################################
    # Write Cell Barcodes, UMIs, and VTs #
#####################################################

vt_whole_fastq = os.path.join(output_dir, vt_f_name)
cb_whole_fastq = os.path.join(output_dir, cb_f_name)
umi_whole_fastq = os.path.join(output_dir, umi_f_name)

print("Writing VTs to: ", vt_whole_fastq)
print("Writing CBCs to: ", cb_whole_fastq)
print("Writing UMIs to: ", umi_whole_fastq)

# files exist remove them
if os.path.exists(vt_whole_fastq):
    os.remove(vt_whole_fastq)

if os.path.exists(cb_whole_fastq):
    os.remove(cb_whole_fastq)

if os.path.exists(umi_whole_fastq):
    os.remove(umi_whole_fastq)

with gzip.open(vt_whole_fastq, 'wt') as f:
    for vt in vts:
        f.write(vt + "\n")

with gzip.open(cb_whole_fastq, 'wt') as f:
    for cb in cellbcs:
        f.write(cb + "\n")

with gzip.open(umi_whole_fastq, 'wt') as f:
    for umi in umis:
        f.write(umi + "\n")


summary_json_f_name = "barcode_extraction_summary.json"
summary_json_f_path = os.path.join(output_dir, summary_json_f_name)

print("Writing summary json to: ", summary_json_f_path)

with open(summary_json_f_path, 'w') as f:
    json.dump(filter_stats, f)

