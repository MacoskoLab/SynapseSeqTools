from datetime import date
today = date.today().strftime("%Y%m%d")

import csv
import gzip

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import scipy.io
import pandas as pd

import os
import glob

from collections import Counter, defaultdict

from shapely.geometry import Point
from shapely.geometry.polygon import Polygon

from itertools import compress
from itertools import combinations

import click

def plot_slideseq_umis(
    slideseq_library, bead_color, bkgnd_color, plot_title, plot_size=(14,10), pct=95
):  
    """
    beads - list of bead barcodes
    bead_xy_d - dictionary from barcode to x,y coordinates
    bead_color - base color for all Slide-seq beads on given puck
    query_gene - name of gene being plotted
    gene_cmap - cmap desired for the gene being plotted
    plot_title - title for the figure
    pct - max percentile for the color scale on the beads, tweak for visualization
    """
        
    fig, ax = plt.subplots(1, 1, figsize=plot_size)
    fig.patch.set_facecolor("white")
    
    umis_data=np.flip(np.asarray(slideseq_library['matrix'].sum(axis=0))[0])
    norm = matplotlib.colors.Normalize(0, np.percentile(umis_data, pct), clip=True)
    
    base = ax.scatter(
        [slideseq_library['bead_xy_d'][b][0] for b in reversed(slideseq_library['barcodes'])],
        [slideseq_library['bead_xy_d'][b][1] for b in reversed(slideseq_library['barcodes'])],
        c=umis_data,
        s=8,
        cmap="viridis_r",
        norm=norm,
    )
    base.set_rasterized(True)

    ax.set_facecolor(bkgnd_color)
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.axis("equal")
    fig.colorbar(base, ax=ax)

    ax.set_title(plot_title)
    
    plt.savefig(plot_title)
    
def plot_slideseq_genes(
    slideseq_library, bead_color, bkgnd_color, query_gene_list, gene_cmap_list, plot_title, plot_size=(14,10), pct=95
):  
    """
    beads - list of bead barcodes
    bead_xy_d - dictionary from barcode to x,y coordinates
    bead_color - base color for all Slide-seq beads on given puck
    query_gene - name of gene being plotted
    gene_cmap - cmap desired for the gene being plotted
    plot_title - title for the figure
    pct - max percentile for the color scale on the beads, tweak for visualization
    """
    if len(query_gene_list) != len(gene_cmap_list):
        sys.exit('Gene list length does not match cmap list length')
        
    fig, ax = plt.subplots(1, 1, figsize=plot_size)
    fig.patch.set_facecolor("white")

    base = ax.scatter(
        [slideseq_library['bead_xy_d'][b][0] for b in slideseq_library['barcodes']],
        [slideseq_library['bead_xy_d'][b][1] for b in slideseq_library['barcodes']],
        color=bead_color,
        s=5,
        alpha=0.5,
    )
    
    gene_plots={}
    for i_gene, i_cmap in zip(query_gene_list,gene_cmap_list):
        gene_plots[i_gene]=plot_slideseq_singleGene(slideseq_library, ax, i_gene, i_cmap, pct=95)
    
    ax.set_facecolor(bkgnd_color)
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.axis("equal")
    
    for i in gene_plots.keys():
        cbar=fig.colorbar(gene_plots[i], ax=ax)
        cbar.ax.set_xlabel(i, rotation=0)
        
    ax.set_title(plot_title)
    
    plt.savefig(plot_title)

def plot_slideseq_singleGene(
    slideseq_library, ax, query_gene, query_cmap, pct=95
):
    """
    beads - list of bead barcodes
    bead_xy_d - dictionary from barcode to x,y coordinates
    bead_color - dictionary from barcode to color data per bead (e.g. nUMIs)
    matched_tags - tags that matched to snSeq
    title - title for the figure
    pct - max percentile for the color scale on the beads, tweak for visualization
    """
    
    query_gene_idx=slideseq_library['features_idx_dict'][query_gene]
    query_gene_idx_idx=np.where(slideseq_library['matrix'].nonzero()[0]==query_gene_idx)[0]
    query_gene_data=slideseq_library['matrix'].data[query_gene_idx_idx]
    query_gene_idx_idx_idx=slideseq_library['matrix'].nonzero()[1][query_gene_idx_idx]
    query_gene_barcodes=[slideseq_library['barcodes_idx_dict'][i] for i in query_gene_idx_idx_idx]
    
    # using the list of beads and dictionaries to get 
    norm = matplotlib.colors.Normalize(0, np.percentile(query_gene_data, pct), clip=True)
    
    curr_plot = ax.scatter(
        [slideseq_library['bead_xy_d'][b][0] for b in query_gene_barcodes],
        [slideseq_library['bead_xy_d'][b][1] for b in query_gene_barcodes],
        c=query_gene_data,
        s=8,
        cmap=query_cmap,
        norm=norm,
    )
    curr_plot.set_rasterized(True)
    
    return curr_plot

def cmap_maker(
    color_list, cmap_name
):
    colors=[matplotlib.colors.to_rgba(i) for i in color_list]
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list(cmap_name, colors, N=1000)
    
    return cmap

def spatial_and_matched(
    beads, bead_xy_d, bead_color, beads_matched,bead_color2, title, plot_title, pct=95
):
    """
    beads - list of bead barcodes
    bead_xy_d - dictionary from barcode to x,y coordinates
    bead_color - dictionary from barcode to color data per bead (e.g. nUMIs)
    matched_tags - tags that matched to snSeq
    title - title for the figure
    pct - max percentile for the color scale on the beads, tweak for visualization
    """
    fig, ax = plt.subplots(1, 1, figsize=(8, 8))
    fig.patch.set_facecolor("white")

    # using the list of beads and dictionaries to get 
    bead_color = np.array([bead_color[b] for b in beads])

    # version of 'Blues' colormap that is pure white at the bottom
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
        "BluesW",
        [(1.0, 1.0, 1.0), (0.0314, 0.188, 0.450)]
    )
    norm = matplotlib.colors.Normalize(0, np.percentile(bead_color, pct), clip=True)

    c = ax.scatter(
        [bead_xy_d[b][0] for b in beads],
        [bead_xy_d[b][1] for b in beads],
        c=bead_color,
        s=0.5,
        alpha=0.5,
        cmap=cmap,
        norm=norm,
    )
    c.set_rasterized(True)
    
    d = ax.scatter(
        [bead_xy_d[b][0] for b in beads_matched],
        [bead_xy_d[b][1] for b in beads_matched],
        c=bead_color2,
        s=12,
        alpha=0.2,
        cmap=cmap,
        norm=norm,
    )
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.axis("equal")
    ax.set_title(title)
    fig.colorbar(c, ax=ax)
    
    plt.savefig(plot_title)    
    
def get_gene_vals(gene):
    gene_idx=dge_features[dge_features==gene].index[0]
    gene_val=pd.Series(dge[gene_idx,:].toarray().flatten())
    
    return gene_val

def dict_to_csr(input_dict, input_cell, input_tags):
    m = scipy.sparse.dok_matrix((len(input_cell), len(input_tags)), dtype=np.int32)
    b2i = {b: i for i, b in enumerate(input_cell)}
    t2j = {t: j for j, t in enumerate(input_tags)}

    for b, bd in input_dict.items():
        for t, v in bd.items():
            m[b2i[b], t2j[t]] = v

    m = m.tocsr()

    return m

def whitelist_filter(input_csr, whitelist, input_cell, input_tags):
    
    temp_data=input_csr.data
    temp_cell=input_csr.nonzero()[0]
    temp_tags=input_csr.nonzero()[1]
    
    whitelist_data=temp_data[whitelist]
    whitelist_cell=temp_cell[whitelist]
    whitelist_tags=temp_tags[whitelist]
    
    whitelist_originalIdx={'data':whitelist_data, 'cell':whitelist_cell, 'tags':whitelist_tags}

    cells_by_tag=defaultdict(lambda: defaultdict(int))
    
    for idx1, idx2, val in zip(whitelist_cell, whitelist_tags, whitelist_data):
        cells_by_tag[input_cell[idx1]][input_tags[idx2]]=val
        
    #cells_by_tag = {
    #        bead_bc: {dtag: len(v) for dtag, v in dtags_per_bead.items()}
    #        for bead_bc, dtags_per_bead in cells_by_tag.items()
    #    }
    
    output_cell = sorted(cells_by_tag)
    output_tags = sorted({t for b in output_cell for t in cells_by_tag[b]})
    output_m = dict_to_csr(cells_by_tag, output_cell, output_tags)
    
    whitelist_newIdx={'m':output_m, 'cell':output_cell, 'tags':output_tags}
    
    return whitelist_originalIdx, whitelist_newIdx

def initial_h_set(barcodes):
    base_d = {"A": 1, "C": 2, "G": 3, "T": 4, "N": 5}
    return [int(''.join(map(str,tuple(base_d[c] for c in barcode)))) for barcode in barcodes]
    
def intersection_idx(base_list,query_list):
    int_idx=np.in1d(initial_h_set(base_list), initial_h_set(query_list))
    return int_idx

def u_to_m(input_u, input_cells, input_tags):
    
    cells_bc = []
    cells_umi = []
    for i in input_cells:
        cells_bc.append(i[:14])
        cells_umi.append(i[14:])
    
    cells_by_tag = defaultdict(lambda: defaultdict(set))
    for idx1, idx2 in zip(input_u.nonzero()[0], input_u.nonzero()[1]):
        cells_by_tag[cells_bc[idx1]][input_tags[idx2]].add(cells_umi[idx1])
        
    cells_by_tag = {
            bead_bc: {dtag: len(v) for dtag, v in dtags_per_bead.items()}
            for bead_bc, dtags_per_bead in cells_by_tag.items()
        }
    
    output_cell = sorted(cells_by_tag)
    output_tags = sorted({t for b in output_cell for t in cells_by_tag[b]})
    output_m = dict_to_csr(cells_by_tag, output_cell, output_tags)
    
    output_all={'m':output_m, 'cell':output_cell, 'tags':output_tags}
    
    return output_all

def polygon_caller(call_path):
    f = open(call_path, "r")
    polygon_raw=[line.strip() for line in f]
    polygon_idx=[(float(i.split(" ")[0]), float(i.split(" ")[1])) for i in polygon_raw]
    polygon_out=Polygon(polygon_idx)
    
    return polygon_out

def polygon_plotter(slideseq_dialout,rmFilt_newIdx,polygon_in,polygon_name,output_name):
    fig, ax = plt.subplots(1, 1, figsize=(8, 8))
    fig.patch.set_facecolor("white")
        
    d = ax.scatter(
        [slideseq_dialout['bead_xy_d'][b][0] for b in rmFilt_newIdx['cell']],
        [slideseq_dialout['bead_xy_d'][b][1] for b in rmFilt_newIdx['cell']],
        c='magenta',
        s=12,
        alpha=0.2,
    )
    ax.axis("equal")

    plt.plot(*polygon_in.exterior.xy,color='blue')  
    plt.title(f"{polygon_name}")
    plt.savefig(output_name)

@click.command()
@click.argument("input_path", type=click.Path(exists=True))
@click.option(
    "--input-sample",
    help="temp",
    required=True,
)
@click.option(
    "--output-folder",
    default="SlideDO_Filter",
    help="temp",
)
@click.option(
    "--polygon-folder",
    help="polygon_maker",
    required=True,
)
@click.option(
    "--polygon-list",
    default="dLGN",
    help="temp",
)
def main(
    input_path,
    input_sample,
    output_folder,
    polygon_folder,
    polygon_list,
):
    polygon_list=polygon_list.split(',')

    curr_path=os.path.join(input_path,input_sample)
    output_path=os.path.join(curr_path,output_folder)
    polygon_path = polygon_folder
    #os.path.join(curr_path,polygon_folder)

    try:
        os.mkdir(output_path)
    except:
        print(f"Output folder exists: {output_path}")

    slideseq_dialout={}
    with gzip.open(os.path.join(curr_path,"raw_umi_matrix.mtx.gz"), "rb") as fh:
        slideseq_dialout['m'] = scipy.io.mmread(fh).tocsr()
    with gzip.open(os.path.join(curr_path,"barcode_matching.txt.gz"), "rt") as fh:
        rows = csv.reader(fh, delimiter="\t")
        slideseq_dialout['bead_xy_d'] = {r[1]: np.array((float(r[2]), float(r[3]))) for r in rows}
    with gzip.open(os.path.join(curr_path,"beads.txt.gz"), "rt") as fh:
        slideseq_dialout['beads'] = [line.strip() for line in fh]
    with gzip.open(os.path.join(curr_path,"raw_tags.txt.gz"), "rt") as fh:
        slideseq_dialout['tags'] = [line.strip() for line in fh]
    with gzip.open(os.path.join(curr_path,"raw_umi_matrix_umi.mtx.gz"), "rb") as fh:
        slideseq_dialout['u'] = scipy.io.mmread(fh).tocsr()
    with gzip.open(os.path.join(curr_path,"beads_umi.txt.gz"), "rt") as fh:
        slideseq_dialout['beads_umi'] = [line.strip() for line in fh]
    slideseq_dialout['bead_numis'] = dict(zip(slideseq_dialout['beads'], np.asarray(slideseq_dialout['m'].sum(1)).flatten()))

    cells_per_tag=Counter(slideseq_dialout['m'].nonzero()[1])
    tmp_idx=np.where(np.asarray(list(cells_per_tag.values()))>1)[0]
    rep_tags_idx=np.asarray(list(cells_per_tag.keys()))[tmp_idx]
    rep_tags = [ slideseq_dialout['tags'][i] for i in rep_tags_idx ]

    spatial_and_matched(
        slideseq_dialout['beads'],
        slideseq_dialout['bead_xy_d'],
        slideseq_dialout['bead_numis'],
        [],
        'magenta',
        f"{len(slideseq_dialout['bead_numis'])} Beads / {sum(np.asarray(list(slideseq_dialout['bead_numis'].values())))} Tags",
        os.path.join(output_path,today+"_1_umiPlot"+'.pdf')
        )

    top5_count=int(np.round(len(slideseq_dialout['u'].data)*0.15))
    rThresh=slideseq_dialout['u'].data[np.argsort(slideseq_dialout['u'].data)[-top5_count:-top5_count+1][0]]
    print(f"Read threshold placed at: {rThresh}")

    rWhitelist=np.arange(0,len(slideseq_dialout['u'].data))[slideseq_dialout['u'].data>=rThresh]
    u_rFilt_ogIdx, u_rFilt_newIdx = whitelist_filter(slideseq_dialout['u'],rWhitelist,slideseq_dialout['beads_umi'],slideseq_dialout['tags'])
    rFilt_newIdx = u_to_m(u_rFilt_newIdx['m'], u_rFilt_newIdx['cell'], u_rFilt_newIdx['tags'])

    ax=spatial_and_matched(
        slideseq_dialout['beads'],
        slideseq_dialout['bead_xy_d'],
        slideseq_dialout['bead_numis'],
        rFilt_newIdx['cell'],
        'magenta',
        f"{len(rFilt_newIdx['cell'])} Beads / {len(rFilt_newIdx['tags'])} Tags",
        os.path.join(output_path,today+"_2_rFilt_slideseq"+'.pdf')
        )

    rFilt_dict=dict( (k,i) for i,k in enumerate(rFilt_newIdx['tags']) )
    mWhitelist_idx=[ rFilt_dict[i] for i in set(rFilt_newIdx['tags']).difference(set(rep_tags)) ]
    mWhitelist=np.where(np.in1d(rFilt_newIdx['m'].nonzero()[1],mWhitelist_idx))[0]

    rmFilt_ogIdx, rmFilt_newIdx = whitelist_filter(rFilt_newIdx['m'],mWhitelist,rFilt_newIdx['cell'],rFilt_newIdx['tags'])

    ax=spatial_and_matched(
        slideseq_dialout['beads'],
        slideseq_dialout['bead_xy_d'],
        slideseq_dialout['bead_numis'],
        rmFilt_newIdx['cell'],
        'magenta',
        f"{len(rmFilt_newIdx['cell'])} Beads / {len(rmFilt_newIdx['tags'])} Tags",
        os.path.join(output_path,today+"_3_rmFilt_slideseq"+'.pdf')
        )

    tmpDict = dict((k,i) for i,k in enumerate(rmFilt_newIdx['cell']))
    polygon_dict={}
    for i in polygon_list:
        try:
            print(i)
            temp_dict={}
            print(polygon_path)
            call_path=os.path.join(polygon_path,i)
            print(call_path)
            temp_dict['polygon'] = polygon_caller(call_path)
            print("ok")

            polygon_plotter(slideseq_dialout,rmFilt_newIdx,temp_dict['polygon'],i,os.path.join(output_path,today+"_4_polygon_"+i+'.pdf'))
        
            temp_dict['polygon_idx']=[temp_dict['polygon'].contains(Point(i)) for i in slideseq_dialout['bead_xy_d'].values()]
            temp_dict['polygon_beads']=np.asarray(list(slideseq_dialout['bead_xy_d'].keys()))[temp_dict['polygon_idx']]
            temp_dict['polygon_beads_filtered']=list(set(rmFilt_newIdx['cell']).intersection(set(temp_dict['polygon_beads'])))
            
            print("ok")

            tmpList = np.asarray([ tmpDict[i] for i in temp_dict['polygon_beads_filtered'] ])
            tmpIdx = np.where(np.in1d(np.asarray(rmFilt_newIdx['m'].nonzero()[0]),tmpList))[0]
            print("ok")
            tmpIdx2 = np.asarray(rmFilt_newIdx['m'].nonzero()[1])[tmpIdx]
            temp_dict['polygon_tags_filtered'] = list({ rmFilt_newIdx['tags'][i] for i in tmpIdx2 })
            print("ok")
            
            tags_per_bead_count=list(Counter(np.asarray(rmFilt_newIdx['m'].nonzero()[0])[tmpIdx]).values())
            print("ok")
            plt.hist(tags_per_bead_count,bins=range(np.max(tags_per_bead_count)+2))
            plt.yscale('log')
            plt.title(f'Identifiers/Beads - {i}')
            plt.xlabel('# of Identifiers')
            plt.ylabel('# of Beads (log scale)')
            plt.savefig(os.path.join(output_path,today+"_5_Identifiers_Per_Bead_"+i+'.pdf'))
            print("ok")
            with gzip.open(os.path.join(output_path,f"{today}_{i}_beads.txt.gz"), "wt") as out:
                print("\n".join(temp_dict['polygon_beads_filtered']), file=out)
            print("ok")
            with gzip.open(os.path.join(output_path,f"{today}_{i}_tags.txt.gz"), "wt") as out:
                print("\n".join(temp_dict['polygon_tags_filtered']), file=out)
            print("ok")

            polygon_dict[i]=temp_dict
        except:
            print(f"Polygon {i} does not exist")
            continue
    print("ok2")
    all_labeled=[]
    for i in polygon_dict.keys():
        all_labeled=all_labeled+polygon_dict[i]['polygon_beads_filtered']
    print("ok2")
    misc_beads_filtered=list(set(rmFilt_newIdx['cell']).difference(set(all_labeled)))
    print("ok2")
    tmpList = np.asarray([ tmpDict[i] for i in misc_beads_filtered ])
    tmpIdx = np.where(np.in1d(np.asarray(rmFilt_newIdx['m'].nonzero()[0]),tmpList))[0]
    tmpIdx2 = np.asarray(rmFilt_newIdx['m'].nonzero()[1])[tmpIdx]
    misc_tags_filtered = list({ rmFilt_newIdx['tags'][i] for i in tmpIdx2 })

    with gzip.open(os.path.join(output_path,f"{today}_Misc_beads.txt.gz"), "wt") as out:
        print("\n".join(misc_beads_filtered), file=out)
    with gzip.open(os.path.join(output_path,f"{today}_Misc_tags.txt.gz"), "wt") as out:
        print("\n".join(misc_tags_filtered), file=out)

if __name__ == "__main__":
    main()

