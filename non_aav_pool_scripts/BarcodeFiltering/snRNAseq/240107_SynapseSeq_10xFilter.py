from datetime import date
today = date.today().strftime("%Y%m%d")

import csv
import gzip
import itertools
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import scipy.io

import click

import os
import glob

import pandas as pd 
import seaborn as sb
from collections import Counter, defaultdict

from itertools import compress
from itertools import combinations

import pickle 

def write_new_matrix(data, cell, tags):
    ce = list(set(cell))
    ta = list(set(tags))
    m = scipy.sparse.dok_matrix((len(ce), len(ta)), dtype=np.int32)
    b2i = {b: i for i, b in enumerate(ce)}
    t2j = {t: j for j, t in enumerate(ta)}

    for n, i in enumerate(data):
        m[b2i[cell[n]], t2j[tags[n]]] = i

    m = m.tocsr()

    return m, ce, ta

def write_matrix(upb, cells, tags):
    m = scipy.sparse.dok_matrix((len(cells), len(tags)), dtype=np.int32)
    b2i = {b: i for i, b in enumerate(cells)}
    t2j = {t: j for j, t in enumerate(tags)}

    for b, bd in upb.items():
        for t, v in bd.items():
            m[b2i[b], t2j[t]] = v

    m = m.tocsr()

    return m
def write_matrix_1(upb, cells, tags):
    m = scipy.sparse.dok_matrix((len(cells), len(tags)), dtype=np.int32)
    b2i = {b: i for i, b in enumerate(cells)}
    t2j = {t: j for j, t in enumerate(tags)}

    for b, bd in upb.items():
        for t, v in bd.items():
            m[b2i[b], t2j[t]] = 1

    m = m.tocsr()

    return m
def write_matrix_lists(data, cells, tags, cells_list, tags_list):
    m = scipy.sparse.dok_matrix((len(cells_list), len(tags_list)), dtype=np.int32)

    for i,j,k in zip(data,cells,tags):
        m[j,k] = i

    m = m.tocsr()

    return m
def initial_h_set(barcodes):
    base_d = {"A": 1, "C": 2, "G": 3, "T": 4, "N": 5}

    return [int(''.join(map(str,tuple(base_d[c] for c in barcode)))) for barcode in barcodes]

def sub_group_10x(list_group_name,list_bcs,query_group):
    t_start=list_group_name.index(query_group)
    t_last=len(list_group_name)-list_group_name[::-1].index(query_group)
    subbed_group=list_bcs[t_start:t_last]
    
    return subbed_group
def intersection_check(cellnum1, query1, cellnum2, query2):
    int_val=len(set(cellnum1[list(set(query1))]).intersection(set(cellnum2[list(set(query2))])))
    return int_val
def intersection_check_3(cellnum1, query1, cellnum2, query2,cellnum3, query3):
    int_val=len(set(cellnum1[list(set(query1))]).intersection(set(cellnum2[list(set(query2))])).intersection(set(cellnum3[list(set(query3))])))
    return int_val
def intersection_check_slideseq_snseq(slideseq_tags, cellnum2, query2):
    a=set(slideseq_tags).intersection(set(cellnum2[list(set(query2))]))
    int_val=len(a)
    return int_val
def intersection_check_slideseq_snseq_val(slideseq_tags, cellnum2, query2):
    a=set(slideseq_tags).intersection(set(cellnum2[list(set(query2))]))
    return list(a)
def read_filt_process(u_in, u_idx, thresh, bc_in, tags_2):
    
    temp_data=u_in.data[u_idx]
    temp_cell=u_in.nonzero()[0][u_idx]
    temp_tags=u_in.nonzero()[1][u_idx]
    
    u2_whitelist_data_filt=temp_data[temp_data >= thresh]
    u2_whitelist_cell_filt=temp_cell[temp_data >= thresh]
    u2_whitelist_tags_filt=temp_tags[temp_data >= thresh]

    cells_by_tag_2=defaultdict(lambda: defaultdict(set))
    
    for idx1, idx2,val in zip(u2_whitelist_cell_filt,u2_whitelist_tags_filt,u2_whitelist_data_filt):
        cells_by_tag_2[bc_in[idx1]][tags_2[idx2]].add(val)
    cells_by_tag_2 = {
            bead_bc: {dtag: len(v) for dtag, v in dtags_per_bead.items()}
            for bead_bc, dtags_per_bead in cells_by_tag_2.items()
        }
    
    cells_new_2 = sorted(cells_by_tag_2)
    tags_new_2 = sorted({t for b in cells_new_2 for t in cells_by_tag_2[b]})
    
    return u2_whitelist_data_filt, u2_whitelist_cell_filt, u2_whitelist_tags_filt, cells_by_tag_2, cells_new_2, tags_new_2

def load_10x_tags(cur_nm, parent_dir):
    do_vars={}
    input_dir = parent_dir+"/{}/".format(cur_nm)
    
    with gzip.open(input_dir+"raw_umi_matrix.mtx.gz", "rb") as fh:
        m = scipy.io.mmread(fh).tocsr()
    with gzip.open(input_dir+"raw_tags.txt.gz", "rt") as fh:
        tags = [line.strip() for line in fh]
    with gzip.open(input_dir+"cells.txt.gz", "rt") as fh:
        cells = [line.strip() for line in fh]
        #cells = [ci[:-] for ci in cells]
    with gzip.open(input_dir+"raw_umi_matrix_umi.mtx.gz", "rb") as fh:
        u = scipy.io.mmread(fh).tocsr()
    with gzip.open(input_dir+"cells_umi.txt.gz", "rt") as fh:
        cells_u = [line.strip() for line in fh]
    do_vars["m"]=m
    do_vars["tags"]=tags
    do_vars["cells"]=cells
    do_vars["u"]=u
    do_vars["cells_u"]=cells_u
    
    return do_vars

def readFilt_10xStrat(library):
    cell_argmax=[i[0] for i in library['u'].argmax(axis=1).tolist()]
    
    new_rows=[i for i in range(library['u'].shape[0])]
    new_cols=cell_argmax
    new_data=[library['u'][i,j] for i,j in zip(new_rows,new_cols)]
    
    return new_rows, new_cols, new_data

def uTOm_readFilter(u_in, thresh, bc_in, tags_2):
    
    out_dict={}
    
    temp_data=u_in.data
    temp_cell=u_in.nonzero()[0]
    temp_tags=u_in.nonzero()[1]


    u2_whitelist_data_filt=temp_data[temp_data > thresh]
    u2_whitelist_cell_filt=temp_cell[temp_data > thresh]
    u2_whitelist_tags_filt=temp_tags[temp_data > thresh]

    cells_by_tag_2=defaultdict(lambda: defaultdict(set))
    
    for idx1, idx2,val in zip(u2_whitelist_cell_filt,u2_whitelist_tags_filt,u2_whitelist_data_filt):
        cells_by_tag_2[bc_in[idx1][:16]][tags_2[idx2]].add(bc_in[idx1][16:])
    
    cells_by_tag_2 = {
        bead_bc: {tag: len(v) for tag, v in tags_per_bead.items()}
        for bead_bc, tags_per_bead in cells_by_tag_2.items()
    }
    cells_new_2 = sorted(cells_by_tag_2)
    tags_new_2 = sorted({t for b in cells_new_2 for t in cells_by_tag_2[b]})
    
    out_dict["u_data_readFilt"]=u2_whitelist_data_filt
    out_dict["u_cell_readFilt"]=u2_whitelist_cell_filt
    out_dict["u_tags_readFilt"]=u2_whitelist_tags_filt
    out_dict["uRF"]=write_matrix(cells_by_tag_2, cells_new_2, tags_new_2)
    out_dict["cell_uRF"]=cells_new_2
    out_dict["tags_uRF"]=tags_new_2
    
    return out_dict

@click.command()
@click.argument("input-csv", type=click.Path(exists=True))

def main(
    input_csv,
):
    """
    """

    ################
    # Loading data #
    ################
    with open(input_csv, "r") as file:
        # Create a CSV reader object
        reader = csv.DictReader(file)

        # Iterate over each row in the CSV file

        i=1
        for row in reader:
            print(f"Reading input CSV")
            print(row.keys())

            parent_dir = row['parent_dir']
            output_dir = row['output_dir']+"/"+today
            target_datasets = row['target_datasets'].split(';')
            custom_prefix = int(row['custom_prefix'])
            cluster_identity_dir = row['cluster_identity_dir']
            custom_cluster = int(row['custom_cluster'])
            perc_cutoff_list = [float(i) for i in row['perc_cutoff_list'].split(';')]

            break
    try:
        os.mkdir(output_dir)
    except:
        print('Output path already exists - continuing')

    print("Processing data")
    if custom_prefix > -1:
        name_prefix = [i.split('_')[custom_prefix] for i in target_datasets]
    else:
        name_prefix = ["S"+i.split('_')[1] for i in target_datasets]
    prefix_2_dataset = {i:j for i,j in zip(name_prefix, target_datasets)}

    print(name_prefix)

    all_libraries_raw = {i:load_10x_tags(i, parent_dir) for i in target_datasets}

    check_cellLen=len(all_libraries_raw[target_datasets[0]]['cells'][0])
    if check_cellLen != 16:
        sys.exit(f'Cell length is incorrect (Expected: 16, Observed: {check_cellLen}')
    else:
        print(f"{all_libraries_raw[target_datasets[0]]['cells'][0]} - Length: {check_cellLen}")

    for i in target_datasets:
        all_libraries_raw[i]["cells_INT"]=np.asarray(initial_h_set(all_libraries_raw[i]["cells"]))
        all_libraries_raw[i]["m_maxTags"]=np.asarray(all_libraries_raw[i]["m"].argmax(axis=0))[0]
        all_libraries_raw[i]["cells_u_bc"]=[j[:16] for j in all_libraries_raw[i]["cells_u"]]
        all_libraries_raw[i]["cells_u_umi"]=[j[16:] for j in all_libraries_raw[i]["cells_u"]]
        all_libraries_raw[i]["cells_INT"]=np.asarray(initial_h_set(all_libraries_raw[i]["cells"]))
        all_libraries_raw[i]["cells_u_bc_INT"]=np.asarray(initial_h_set(all_libraries_raw[i]["cells_u_bc"]))    

    #################
    # 10x Filtering #
    #################

    print('Filtering SSI - 10x Strategy')

    for i in target_datasets:
        new_rows, new_cols, new_data = readFilt_10xStrat(all_libraries_raw[i])
        all_libraries_raw[i]['u_10xFilt'] = write_matrix_lists(new_data, new_rows, new_cols, all_libraries_raw[i]['cells_u'], all_libraries_raw[i]['tags'])

    ######################
    # Cluster Separation #
    ######################

    print("Separating data by transcriptomic cell clusters") 
    cluster_identity_paths=glob.glob(cluster_identity_dir+"/*.csv")
    if custom_cluster > -1:
        cluster_names = [i.split("/")[-1].split(".")[0].split("_")[custom_cluster] for i in cluster_identity_paths]
    else:
        cluster_names = [i.split("/")[-1].split(".")[0].replace("_","") for i in cluster_identity_paths]

    # Opening the files for each cluster
    cluster_by_barcodes={}
    temp_set = set()
    for i,j in zip(cluster_identity_paths, cluster_names):
        with open(i) as fh:
            csvfh = csv.reader(fh)
            next(csvfh)
            cluster_by_barcodes[j]=[entry.split('_')[1].split("-")[0] for row,entry in csvfh]
        with open(i) as fh:
            csvfh = csv.reader(fh)
            next(csvfh)
            cluster_by_barcodes[j+"_Dataset"]=[entry.split('_')[0] for row,entry in csvfh]
        temp_set = temp_set.union(set(cluster_by_barcodes[j+"_Dataset"]))

    #name_prefix = sorted(list(temp_set))
    #prefix_2_dataset = {i:j for i,j in zip(name_prefix, target_datasets)}


    # Separating clusters by their originating dataset
    cluster_by_barcodes_by_dataset={}
    for i in cluster_names:
        for grpcode in name_prefix:
            try:
                cluster_by_barcodes_by_dataset[grpcode+"_"+i]=sub_group_10x(cluster_by_barcodes[i+"_Dataset"], cluster_by_barcodes[i], grpcode)
            except:
                continue

    # Getting the index for the cell barcodes / cell+UMI barcodes for each cluster
    cluster_by_barcodes_by_dataset_idx={}
    cluster_by_barcodes_by_dataset_u_idx={}
    for i in cluster_by_barcodes_by_dataset.keys():
        j=i.split('_')[0] # The group code
        k=i.split('_')[1] # The actual cluster
        
        cluster_by_barcodes_by_dataset_idx[i+"_idx"] = np.in1d(all_libraries_raw[prefix_2_dataset[j]]["cells_INT"], initial_h_set(cluster_by_barcodes_by_dataset[i]))
        cluster_by_barcodes_by_dataset_u_idx[i+"_idx"] = np.in1d(all_libraries_raw[prefix_2_dataset[j]]["cells_u_bc_INT"], initial_h_set(cluster_by_barcodes_by_dataset[i]))
    
    # Converting clusters into essentially separate datasets
    all_libraries_byCluster_u_noFilt={}
    for i in cluster_by_barcodes_by_dataset_u_idx.keys():
        j=i.split('_')[0] # The group code
        
        m_temp=all_libraries_raw[prefix_2_dataset[j]]["u_10xFilt"]
        temp=cluster_by_barcodes_by_dataset_u_idx[i]
        
        all_libraries_byCluster_u_noFilt[i+"_cell"] = m_temp.nonzero()[0][np.in1d(m_temp.nonzero()[0],np.asarray(list(range(0,len(temp))))[temp])]
        all_libraries_byCluster_u_noFilt[i+"_tags"] = m_temp.nonzero()[1][np.in1d(m_temp.nonzero()[0],np.asarray(list(range(0,len(temp))))[temp])]
        all_libraries_byCluster_u_noFilt[i+"_umis"] = m_temp.data[np.in1d(m_temp.nonzero()[0],np.asarray(list(range(0,len(temp))))[temp])]

    ##################
    # Read Filtering #
    ##################

    print('Filtering by reads')
    for perc_cutoff in perc_cutoff_list:
        perc_cutoff_name = int(100*perc_cutoff)
        output_folder = f'{output_dir}/Perc{perc_cutoff_name}/'

        try:
            os.mkdir(output_folder)
        except:
            print('Output path already exists - continuing')

        # Calculating thresholds and plotting
        thresholds={}
        for i in target_datasets:
            percentile_value = int(np.round(len(all_libraries_raw[i]['u_10xFilt'].data)*perc_cutoff))
            threshold_value = all_libraries_raw[i]['u_10xFilt'].data[np.argsort(all_libraries_raw[i]['u_10xFilt'].data)[-percentile_value:-percentile_value+1][0]]

            thresholds[i]=threshold_value

        fig, axs = plt.subplots(len(target_datasets))
        for n,i in enumerate(target_datasets):

            tmp=all_libraries_raw[i]["u_10xFilt"].data

            axs[n].hist(tmp,bins=range(0,300,1))
            axs[n].axvline(x=thresholds[i]+1,color='r')
            axs[n].set_title("{} - Thresh: {}".format(i,thresholds[i]))
            axs[n].set_yscale('log') #, nonposy='clip')

        fig.set_figwidth(8)
        fig.set_figheight(31)
      
        fig.savefig(output_folder+today+"_ReadThresh_10xDO_Histograms.png")

        # Actually filtering and subsetting
        all_libraries_READFILT={i: uTOm_readFilter(all_libraries_raw[i]["u_10xFilt"], thresholds[i], all_libraries_raw[i]["cells_u"], all_libraries_raw[i]["tags"]) for i in target_datasets}

        all_libraries_READFILT_byCluster={}
        for i in cluster_by_barcodes_by_dataset_u_idx.keys():
            j=i.split('_')[0] # the group code
            
            temp = np.in1d(all_libraries_READFILT[prefix_2_dataset[j]]["u_cell_readFilt"],all_libraries_byCluster_u_noFilt[i+"_cell"])

            all_libraries_READFILT_byCluster[i+"_cell"] = all_libraries_READFILT[prefix_2_dataset[j]]["u_cell_readFilt"][temp]
            all_libraries_READFILT_byCluster[i+"_tags"] = all_libraries_READFILT[prefix_2_dataset[j]]["u_tags_readFilt"][temp]
            all_libraries_READFILT_byCluster[i+"_data"] = all_libraries_READFILT[prefix_2_dataset[j]]["u_data_readFilt"][temp]

        # Plotting UMIs/Cell Type boxplots
        counter_list={}
        for i in cluster_by_barcodes_by_dataset_u_idx.keys():
            k=i.split('_')[0]
            
            counter_list[i]=list(Counter([all_libraries_raw[prefix_2_dataset[k]]["cells_u_bc"][j] for j in all_libraries_READFILT_byCluster[i+"_cell"]]).values())

        df = pd.DataFrame(columns = ['Celltypes', 'IDs'])
        for key, value in counter_list.items():
            temp_df = pd.DataFrame({'Celltypes': [key.split('_')[1]]*len(value),
                                       'IDs': [float(i) for i in value]})
            df = df.append(temp_df, ignore_index = True)

        fig, axes = plt.subplots()

        sb.boxplot(x='Celltypes',y='IDs', data=df, ax = axes)
        axes.set_title('Synapse-seq IDs per Cell')

        axes.yaxis.grid(True)
        axes.set_yscale('log') #, nonposy='clip')
        axes.set_ylim([1, 10000])
        #axes.set_xlabel('Scenario')
        #axes.set_ylabel('LMP ($/MWh)')

        fig.set_figwidth(16)
        fig.set_figheight(8)

        plt.show()

        fig.savefig(output_folder+today+"_ids_per_cell.png")

        # Actually saving the data
        all_libraries_READFILT_byCluster_Output={}
        for i in cluster_by_barcodes_by_dataset_u_idx.keys():
            j=i.split('_')[0]
            
            all_libraries_READFILT_byCluster_Output[i+"_cell"]=[all_libraries_raw[prefix_2_dataset[j]]["cells_u_bc"][i] for i in all_libraries_READFILT_byCluster[i+"_cell"]]
            all_libraries_READFILT_byCluster_Output[i+"_tags"]=[all_libraries_raw[prefix_2_dataset[j]]["tags"][i] for i in all_libraries_READFILT_byCluster[i+"_tags"]]
            all_libraries_READFILT_byCluster_Output[i+"_data"]=all_libraries_READFILT_byCluster[i+"_data"]

        with open(output_folder+today+'_readFilter_10xDO_Data.pkl', 'wb') as f:
            pickle.dump(all_libraries_READFILT_byCluster_Output, f)

        for i in cluster_by_barcodes_by_dataset_u_idx.keys():

            m,ce,ta = write_new_matrix(all_libraries_READFILT_byCluster_Output[i+"_data"], all_libraries_READFILT_byCluster_Output[i+"_cell"], all_libraries_READFILT_byCluster_Output[i+"_tags"])

            output_tags = output_folder+today+'_'+i+'_tags.txt.gz'
            with gzip.open(output_tags, "wt") as out:
                print("\n".join(ta), file=out)

            output_cell = output_folder+today+'_'+i+'_cells.txt.gz'
            with gzip.open(output_cell, "wt") as out:
                print("\n".join(ce), file=out)

            output_file = output_folder+today+'_'+i+'_matrix.mtx.gz'
            with gzip.open(output_file, "wb") as out:
                scipy.io.mmwrite(out, m)

if __name__ == "__main__":
    main()
