import gzip
import itertools
from collections import Counter, defaultdict
import pickle

import csv
import os
import numpy as np
import scipy.io
import scipy.sparse
from sklearn.neighbors import BallTree
import sys

import argparse
from matplotlib import pyplot as plt
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
import multiprocessing as mp
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
sys.path.append("/home/jsilverm/06_synapseseq_repo/synapse_seq_pipeline_code/09_dedup_aav_lib/")
import synapse_seq_functions as ssf


def polygon_caller(call_path):
    f = open(call_path, "r")
    polygon_raw=[line.strip() for line in f]
    polygon_idx=[(float(i.split(" ")[0]), float(i.split(" ")[1])) for i in polygon_raw]
    polygon_out=Polygon(polygon_idx)
    
    return polygon_out

def create_region_to_vt_mapping(polygon_folder, polygon_list, read_filt_obj,bead_coord_dict):
    region_to_vt_set = {}
    assert read_filt_obj.keys() == set({"u_mat", "m_mat", "cells", "vts", "cells_umi"})
    for region_name in polygon_list:
        polygon_coords_path = os.path.join(polygon_folder, region_name)
        assert os.path.exists(polygon_coords_path)
        polygon = polygon_caller(polygon_coords_path)
        polygon_idx = [polygon.contains(Point(i)) for i in bead_coord_dict.values()]
        polygon_beads = np.asarray(list(bead_coord_dict.keys()))[polygon_idx]
        polygon_beads_filtered = set(read_filt_obj["cells"]).intersection(set(polygon_beads))
        polygon_beads_filtered_indices = [i for i, bead_bc in enumerate(read_filt_obj["cells"]) if bead_bc in polygon_beads_filtered]
        m_mat = read_filt_obj["m_mat"]
        # for bead filtered inx, take all nonzero vts
        m_mat_cells_in_region = m_mat[polygon_beads_filtered_indices, :]
        # get nonzero col elements
        nonzero_cols = m_mat_cells_in_region.nonzero()[1]
        vts_of_nonzero_cols = np.array(read_filt_obj["vts"])[nonzero_cols]
        vts_observed_set = set(vts_of_nonzero_cols)

        res = {"region_name": region_name, "vts": vts_observed_set, "bead_names": polygon_beads_filtered, "polygon_coords_path": polygon_coords_path}

        region_to_vt_set[region_name] = res 

    # also add a region for all beads
    nonzero_cols = read_filt_obj["m_mat"].nonzero()[1]
    vts_of_nonzero_cols = np.array(read_filt_obj["vts"])[nonzero_cols]
    vts_observed_set = set(vts_of_nonzero_cols)
    region_to_vt_set["all_beads"] = {"region_name": "all_beads", "vts": vts_observed_set, "bead_names": set(read_filt_obj["cells"])}


    return region_to_vt_set

def plot_beads_in_polygon(polygon_coords_path, polygon_name, read_filt_obj,bead_coord_dict, ax=None):

    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(8, 8))
        # fig.patch.set_facecolor("white")

    polygon_in = polygon_caller(polygon_coords_path)

        
    d = ax.scatter(
        [bead_coord_dict[b][0] for b in read_filt_obj["cells"]],
        [bead_coord_dict[b][1] for b in read_filt_obj['cells']],
        c='magenta',
        s=12,
        alpha=0.2,
    )

    ax.plot(*polygon_in.exterior.xy,color='blue')  
    ax.set_title(f"{polygon_name}")


# working_dir="/home/jsilverm/06_synapseseq_repo/data/00_jonah_processing/03_slideseq_barcode_filtering/dLGN_01/"
# bc_processing_dir="vt_correction"
# polygon_folder="polygon_maker"
# polygon_list="l_dLGN;ll_dLGN;m_dLGN;mm_dLGN"
# bead_umi_split_inx=14

# read args from command line

args_parser = argparse.ArgumentParser()
args_parser.add_argument("-w", "--working_dir", help="Working directory")
args_parser.add_argument("-b", "--bc_processing_dir", help="Barcode processing directory")
args_parser.add_argument("-l", "--polygon_list", help="Polygon list")
args_parser.add_argument("-p", "--polygon_folder", default="polygon_maker", help="Polygon folder")
args_parser.add_argument("-c", "--bc_coords_fname", default="barcode_coordinates.txt.gz" ,help="path to coords of slide-seq beads")


args = args_parser.parse_args()

working_dir = args.working_dir
bc_processing_dir = args.bc_processing_dir
polygon_folder = args.polygon_folder
polygon_list = args.polygon_list
# bc_coords_path=args.bc_coords_path

print(f"Working dir: {working_dir}")
print(f"Barcode processing dir: {bc_processing_dir}")
print(f"Polygon folder: {polygon_folder}")
print(f"Polygon list: {polygon_list}")

sys.stderr.flush()
sys.stdout.flush()

delimiter=";"
polygon_list = polygon_list.split(delimiter)


# polygon_folder_full = os.path.join(working_dir, polygon_folder)
bc_coords_path = os.path.join(working_dir, "barcode_coordinates.txt.gz")
assert os.path.exists(bc_coords_path)

## Hardcoded for now
bead_umi_split_inx = 14
reads_umi_thresholds = [0, 10, 50, 100, 500]

# Parse slideseq beads array
bead_coord_dict = {}
with gzip.open(bc_coords_path, "rt") as fh:
    for line in fh:
        bead, x, y = line.strip().split("\t")
        bead_coord_dict[bead] = (float(x), float(y))

print("Read in bead coordinates")


polygon_folder_full = os.path.join(working_dir, polygon_folder)
bc_processing_folder_full = os.path.join(working_dir, bc_processing_dir)

assert os.path.exists(polygon_folder_full)
assert os.path.exists(bc_processing_folder_full)

# Read in u_mat, filter, and create m_mat
u_mat_obj_path = os.path.join(bc_processing_folder_full, "dedup_vt_and_umi_obj.pkl")
assert os.path.exists(u_mat_obj_path)

print("Reading in dedup object")
with open(u_mat_obj_path, "rb") as fh:
    u_mat_obj = pickle.load(fh)

sys.stderr.flush()
sys.stdout.flush()


# create read thresholded objects at each threshold
thresholded_objs = {}
for thresh in reads_umi_thresholds:
    print(thresh)
    sys.stderr.flush()
    sys.stdout.flush()

    thresholded_obj = ssf.create_u_and_m_mat_from_read_filter(
        u_mat_obj["u"], 
        cellbc_umi_ordering = u_mat_obj["bead_umi"],
        vt_ordering = u_mat_obj["tags"],
        threshold_value =  thresh, 
        cb_split_inx = bead_umi_split_inx)

    thresholded_objs[thresh] = thresholded_obj


region_to_vts_different_filts = {}
for thresh, thresholded_obj in thresholded_objs.items():
    region_to_vts_different_filts[thresh] = create_region_to_vt_mapping(polygon_folder_full, polygon_list, thresholded_obj, bead_coord_dict)


# save region to vts mapping
region_to_vts_out_path = os.path.join(bc_processing_folder_full, "region_to_vts.pkl")
with open(region_to_vts_out_path, "wb") as fh:
    pickle.dump(region_to_vts_different_filts, fh)


# for each thresholded object, create a folder, pickle the object and save it
for thresh, thresholded_obj in thresholded_objs.items():
    out_folder = os.path.join(bc_processing_folder_full, f"reads_{thresh}")
    if not os.path.exists(out_folder):
        os.makedirs(out_folder)

    out_path_pkl = os.path.join(out_folder, "read_filt_dedup_obj.pkl")
    with open(out_path_pkl, "wb") as fh:
        pickle.dump(thresholded_obj, fh)


pdf_out_name = "slide_seq_region_match_read_filt.pdf"
out_pdf_path = os.path.join(bc_processing_folder_full, pdf_out_name)
pp = PdfPages(out_pdf_path)


for thresh, thresholded_obj in thresholded_objs.items():
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5))

    # create umi_tag plot
    ssf.create_umi_per_tag_histogram(thresholded_obj["u_mat"], title_base="Reads > " + str(thresh), ax=ax1)
    ssf.create_cells_per_tag_histogram(thresholded_obj["m_mat"], title_base="Reads > " + str(thresh), ax=ax2)
    ssf.create_tags_per_cell_histogram(thresholded_obj["m_mat"], title_base="Reads > " + str(thresh), ax=ax3)

    pp.savefig(fig)


n_regions = len(polygon_list)
for i in range(n_regions):
    region_name = polygon_list[i]
    region_to_vts = region_to_vts_different_filts[0][region_name]
    polygon_coords_path = region_to_vts["polygon_coords_path"]
    fig, ax1 = plt.subplots(1, 1, figsize=(8, 8))
    plot_beads_in_polygon(polygon_coords_path, region_name, thresholded_objs[0], bead_coord_dict, ax=ax1)
    pp.savefig(fig)



pp.close()