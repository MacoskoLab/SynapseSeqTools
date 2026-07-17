import gzip
import itertools
import logging
from collections import Counter, defaultdict
from pathlib import Path
import pickle

import click
import matplotlib
import matplotlib.colors
import networkx as nx
import numpy as np
import scipy.io
import scipy.sparse
from matplotlib.backends.backend_pdf import PdfPages
from sklearn.neighbors import radius_neighbors_graph

import slideseq.bead_matching
from slideseq.plot import new_ax
from slideseq.util.logger import create_logger

import subprocess
import os

log = logging.getLogger("barcode_matrix")

DEGENERATE_BASE_DICT = {
    "A": {"A"},
    "C": {"C"},
    "G": {"G"},
    "T": {"T"},
    "R": {"A", "G"},
    "Y": {"C", "T"},
    "M": {"A", "C"},
    "K": {"G", "T"},
    "S": {"C", "G"},
    "W": {"A", "T"},
    "H": {"A", "C", "T"},
    "B": {"C", "G", "T"},
    "V": {"A", "C", "G"},
    "D": {"A", "G", "T"},
    "N": {"A", "C", "G", "T"},
}

def match_tags(
    cell_ID, cell_BC, umis, barcode_codes, barcode_cutoff
):
    no_id_struct = 0
    pass_all = 0

    raw_umis_per_bead = defaultdict(lambda: defaultdict(set))
    raw_umis_per_bead_umi = defaultdict(Counter)

    for tag, bead_bc, umi in zip(cell_ID, cell_BC, umis):

        # Skip if Cell ID does not match structure
        if sum(b in s for b, s in zip(tag, barcode_codes)) < barcode_cutoff:
            no_id_struct += 1
            continue

        pass_all += 1

        bead_umi_bc = bead_bc + umi

        raw_umis_per_bead[bead_bc][tag].add(umi)
        raw_umis_per_bead_umi[bead_umi_bc][tag] += 1

    log.debug(f"IDs w/o correct ID structure: {no_id_struct}")
    log.debug(f"Passed all: {pass_all}")

    raw_umis_per_bead = {
        bead_bc: {tag: len(v) for tag, v in tags_per_bead.items()}
        for bead_bc, tags_per_bead in raw_umis_per_bead.items()
    }

    return raw_umis_per_bead, raw_umis_per_bead_umi


def write_matrix(upb, beads, tags, output_file):
    m = scipy.sparse.dok_matrix((len(beads), len(tags)), dtype=np.int32)
    b2i = {b: i for i, b in enumerate(beads)}
    t2j = {t: j for j, t in enumerate(tags)}

    for b, bd in upb.items():
        for t, v in bd.items():
            m[b2i[b], t2j[t]] = v

    m = m.tocsr()

    with gzip.open(output_file, "wb") as out:
        scipy.io.mmwrite(out, m)

    return m

def reverse_h_set(h_barcodes):
    base_d_rev = {0: "A", 1: "C", 2: "G", 3: "T", 4: "N"}

    return ["".join([ base_d_rev[c] for c in h_list ]) for h_list in h_barcodes]

@click.command()
@click.argument("fastq-r1", type=click.Path(exists=True))
@click.argument("fastq-r2", type=click.Path(exists=True))
@click.option(
    "--output-dir",
    help="Path for output",
    required=True,
)
@click.option(
    "--tag-sequence",
    default="WSBWSDWSHWSVWSBWSWSHWSVWSBWSDWSH",
    help="Format for the tag sequence",
)
@click.option(
    "--constant-sequence-wpre",
    default="GATACCGAGCGCTGC",
    help="Constant sequence that should precede the tag",
)
@click.option(
    "--constant-sequence-polya",
    default="TCGAGAGATCTACGGG",
    help="Constant sequence that should follow the tag",
)
@click.option(
    "--tag-mismatch",
    type=int,
    default=3,
    help="Number of mismatches to allow in tag",
    show_default=True,
)
@click.option(
    "--const-mismatch",
    type=int,
    default=3,
    help="Number of mismatches to allow in tag",
    show_default=True,
)
@click.option(
    "--run-merge",
    default=True,
    help="Number of mismatches to allow in tag",
    show_default=True,
)
@click.option(
    "--merge-path",
    default="/broad/macosko/mkim/fake.fastq.gz",
    help="Path for output",
)
@click.option("--debug", is_flag=True, help="Turn on debug logging")
def main(
    fastq_r1,
    fastq_r2,
    output_dir,
    tag_sequence,
    constant_sequence_wpre,
    constant_sequence_polya,
    tag_mismatch,
    const_mismatch,
    run_merge,
    merge_path,
    debug=False,
):
    """
    This script generates some plots for mapping barcoded reads.

    Reads sequences from FASTQ_R1 and FASTQ_R2. Assumes that the first read
    contains a 15bp barcode split across two locations, along with an 8bp UMI.
    The second read is assumed to have TAG_SEQUENCE in bases 20-40.
    """
    create_logger(debug, dryrun=False)

    try:
        os.mkdir(output_dir)
    except:
        print('Path already exists')


    output_dir = Path(output_dir)
    log.info(f"Saving output to {output_dir}")

    command = [f'zcat {fastq_r1} | head -n 2 | tail -n 1 | wc -c']
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    rd1_length = result.stdout.strip()

    rd1_length=int(rd1_length)
    a=rd1_length-1
    log.info(f'Read 1 is {a} nt long')

    if run_merge:
        log.info("Saving only the sequence lines from read 1")
        command = [f'zcat {fastq_r1} | awk \'NR%4==2\' > rd1_stripped.fastq']
        result = subprocess.run(command, cwd=output_dir, shell=True, capture_output=True)

        log.info("Saving only the sequence lines from read 2")
        command = [f'zcat {fastq_r2} | awk \'NR%4==2\' > rd2_stripped.fastq']
        result = subprocess.run(command, cwd=output_dir, shell=True, capture_output=True)

        log.info("Merging Read 1 and Read 2 files")
        command = [f'paste rd1_stripped.fastq rd2_stripped.fastq | gzip > merged.fastq.gz']
        #command = [f'paste -d "" <(zcat {fastq_r1}) <(zcat {fastq_r2}) | awk \'NR%4==2\' | gzip > merged.fastq.gz']
        result = subprocess.run(command, cwd=output_dir, shell=True, capture_output=True)
    
    else:
        log.info("Copying merge file to new folder")
        command = [f'cp {merge_path} merged.fastq.gz']
        result = subprocess.run(command, cwd=output_dir, shell=True, capture_output=True)


    log.info("Fuzzy grep-ing for WPRE constant sequence")
    print(constant_sequence_wpre)
    command = [f'zcat merged.fastq.gz | /broad/macosko/mkim/Scripts/ugrep/bin/ug -Z3 \"{constant_sequence_wpre}\" | gzip > merged_fuzzyMatch.fastq.gz']
    result = subprocess.run(command, cwd=output_dir, shell=True, capture_output=True, text=True)

    log.info("Fuzzy grep-ing for PolyA constant sequence")
    print(constant_sequence_polya)
    command = [f'zcat merged_fuzzyMatch.fastq.gz | /broad/macosko/mkim/Scripts/ugrep/bin/ug -Z3 \"{constant_sequence_polya}\" | gzip > merged_fuzzyMatch2.fastq.gz']
    result = subprocess.run(command, cwd=output_dir, shell=True, capture_output=True, text=True)

    start_id=int(rd1_length)+26
    end_id=int(rd1_length)+57
    log.info("Cutting out the barcodes/identifiers/UMIs")
    command = [f'zcat merged_fuzzyMatch2.fastq.gz | cut -c {start_id}-{end_id} | gzip > merged_fuzzyMatch_notMatch_cellID.fastq.gz; zcat merged_fuzzyMatch2.fastq.gz | cut -c 1-16 | gzip > merged_fuzzyMatch_notMatch_cellBC.fastq.gz; zcat merged_fuzzyMatch2.fastq.gz | cut -c 17-28 | gzip > merged_fuzzyMatch_notMatch_UMI.fastq.gz']
    result = subprocess.run(command, cwd=output_dir, shell=True, capture_output=True, text=True)


    log.info(f"Reading barcodes/identifiers/UMIs")
    cellID_path = output_dir / 'merged_fuzzyMatch_notMatch_cellID.fastq.gz'
    cellBC_path = output_dir / 'merged_fuzzyMatch_notMatch_cellBC.fastq.gz'
    umi_path = output_dir / 'merged_fuzzyMatch_notMatch_UMI.fastq.gz'
    with gzip.open(cellID_path, "rt") as fh:
        cell_ID = [line.strip() for line in itertools.islice(fh, 0, None, 1)]
    with gzip.open(cellBC_path, "rt") as fh:
        cell_BC = [line.strip() for line in itertools.islice(fh, 0, None, 1)]
    with gzip.open(umi_path, "rt") as fh:
        umis = [line.strip() for line in itertools.islice(fh, 0, None, 1)]

    assert len(cell_ID) == len(cell_BC), "read different number of reads"
    assert len(cell_ID) == len(umis), "read different number of reads"

    log.info(f"Total of {len(cell_ID)} reads")

    barcode_codes = [DEGENERATE_BASE_DICT[b] for b in tag_sequence]
    barcode_cutoff = len(tag_sequence) - tag_mismatch
    log.info(f"barcode_cutoff: {barcode_cutoff}")


    raw_umis_per_bead, raw_umis_per_bead_umi = match_tags(
        cell_ID, cell_BC, umis, barcode_codes, barcode_cutoff
    )


    log.info("Writing output files")
    cells = sorted(raw_umis_per_bead)
    cells_umi = sorted(raw_umis_per_bead_umi)
    raw_tags = sorted({t for b in cells for t in raw_umis_per_bead[b]})

    log.info(f"{len(cells)} unique cell barcodes")
    log.info(f"{len(cells_umi)} unique cell barcodes + UMIs")
    log.info(f"{len(raw_tags)} unique cell IDs")

    with gzip.open(output_dir / "cells.txt.gz", "wt") as out:
        print("\n".join(cells), file=out)

    with gzip.open(output_dir / "cells_umi.txt.gz", "wt") as out:
        print("\n".join(cells_umi), file=out)

    with gzip.open(output_dir / "raw_tags.txt.gz", "wt") as out:
        print("\n".join(raw_tags), file=out)

    log.debug("Writing raw umi matrix")
    m = write_matrix(raw_umis_per_bead, cells, raw_tags, output_dir / "raw_umi_matrix.mtx.gz")

    log.debug("Writing raw umi matrix umi")
    write_matrix(raw_umis_per_bead_umi, cells_umi, raw_tags, output_dir / "raw_umi_matrix_umi.mtx.gz")

    log.info("Done!")


if __name__ == "__main__":
    main()


