#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import logging
import os
import re
import shutil
import sys
import time
from multiprocessing import Pool, cpu_count
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO
from python_algorithms.basic.union_find import UF

from common.align_fasta import FastaAligner
from common.prepare_fasta import PrepFasta as pf
from common.prime_fasta import PrimeFastaWriter
from common.prime_fasta_processing import FastaFinalizer

# TODO
# [] make proper help
# [] replace all simple paths by PurePath
# [] спросить Катю про KK17CL172


parser = argparse.ArgumentParser()
parser.add_argument("-i", help="path to RE results (top level)",
                    type=str,
                    metavar="path")
parser.add_argument("-p", help="prefix for fasta filename",
                    type=str,
                    metavar="prefix")
parser.add_argument("out", help="path to results directory")
parser.add_argument("known", help="fasta with known repeats")
parser.add_argument("-l", help="save logfile in output directory",
                    action="store_true", dest="log")
parser.add_argument("-c", help="number of CPU to use",
                    type=int, default=cpu_count(),
                    dest="cpu_number", metavar="CPU")
parser.add_argument("-eo", "--exclude-other",
                    help=("exclude all 'other' contigs and clusters from \
                        analysis (default: True)"),
                    action="store_true", dest="exclude_other")
parser.add_argument("-ir", "--include-ribosomal",
                    action="store_false",
                    help="include rDNA clusters (rank 4) in analysis (default: False)",
                    dest="include_ribosomal")
args = parser.parse_args()

# create folder structure
out_path = Path(args.out)
out_path.mkdir(parents=True, exist_ok=True)

fasta_path = Path(args.out).joinpath("fasta")
fasta_path.mkdir(parents=True, exist_ok=True)

prime_fasta = Path(args.out).joinpath("results", "prime_fasta")
prime_fasta.mkdir(parents=True, exist_ok=True)

final_fasta = Path(args.out).joinpath("results", "final_fasta")
final_fasta.mkdir(parents=True, exist_ok=True)

# logging
if args.log:
    logging.basicConfig(level=logging.DEBUG,
                        filename=args.out + "/REcomp.log",
                        format="\n%(asctime)s - %(funcName)s - %(levelname)s -\n%(message)s\n",
                        filemode="w")
    console = logging.StreamHandler()
    console.setLevel(logging.DEBUG)
    # console.setLevel(logging.INFO)
    formatter = logging.Formatter("\n%(asctime)s - %(funcName)s - %(levelname)s -\n%(message)s\n")
    console.setFormatter(formatter)
    logging.getLogger("").addHandler(console)
else:
    logging.basicConfig(level=logging.INFO,
                        format="\n%(asctime)s - %(funcName)s - %(levelname)s -\n%(message)s\n")

logging.info("""
------------------------------------------------------------------
PIPELINE VERSION               :           2.0.0

AUTHOR                         :      Aleksey Ermolaev
------------------------------------------------------------------
             """)
logging.info(args)
# prepare fasta files with ranks and "others"
work_dirs = {path:prefix for path, prefix in zip(args.i.split(),
                                                 args.p.split())}
# fasta_to_unite = []
if args.exclude_other:
    logging.info("creating fasta containing all contigs, ranks and references")
else:
    logging.info("creating fasta containing all ranks and references")
for path, prefix in work_dirs.items():
    fasta_prep = pf(path, args.known, prefix)
    fasta_prep.create_united_fasta(fasta_path,
                                   exclude_other=args.exclude_other,
                                   include_ribosomal=args.include_ribosomal)
with open(Path(fasta_path).joinpath("fasta.fasta"), "a") as fasta:
    for record in SeqIO.parse(args.known, "fasta"):
        SeqIO.write(record, fasta, "fasta")

records_number = 0
# chunk fasta for parallel
record_iter = SeqIO.parse(open(Path(fasta_path).joinpath("fasta.fasta")), "fasta")
for i, batch in enumerate(fasta_prep.batch_iterator(record_iter, 1000)):
    records_number += len(batch)
    filename = Path(fasta_path).joinpath(f"fasta{i}.fasta")
    with open(filename, "w") as handle:
        count = SeqIO.write(batch, handle, "fasta")
    print(f"Wrote {count} records to {filename.name}")

# prepare connectivity table
fasta_aligner = FastaAligner(fasta_path.joinpath("fasta.fasta"))
files = [path for path in fasta_path.rglob("*.fasta")
         if any(map(str.isdigit, Path(path).stem))]
pool = Pool(processes=args.cpu_number)
result = pool.map(fasta_aligner.align_fasta, files)
blast_table = pd.concat(result)
pool.close()
blast_table = blast_table[blast_table["qseqid"] != blast_table["sseqid"]]
blast_table.sort_values(["qseqid", "sseqid"], 0,
                        inplace=True, ignore_index=True)
blast_table.drop_duplicates(keep="first", inplace=True, ignore_index=True)
con_table = list(blast_table.to_records(index=False))

# prepare data for UF
map_dict = {}
counter = 0
for record in SeqIO.parse(Path(fasta_path).joinpath("fasta.fasta"), "fasta"):
    map_dict[record.id] = counter
    counter += 1
blast_table["qseqid"] = blast_table["qseqid"].map(map_dict)
blast_table["sseqid"] = blast_table["sseqid"].map(map_dict)

# quick union
quick_union = UF(len(map_dict))
for pair in blast_table.itertuples(index=False, name=None):
    quick_union.union(pair[0], pair[1])
uf_repr = [int(i) for i in str(quick_union).split()]
cc_num = set(uf_repr)

# create prime fasta with superclusters
logging.info("prepare primary fasta")
rev_map_dict = {v:k for k,v in map_dict.items()} # invert map dict
prime_fw = PrimeFastaWriter(fasta_path.joinpath("fasta.fasta"),
                            prime_fasta,
                            final_fasta,
                            rev_map_dict,
                            uf_repr)
pool = Pool(processes=args.cpu_number)
pool.map(prime_fw.write_fasta, cc_num)
pool.close()


# process prime fasta into final
fasta = [path for path in prime_fasta.rglob("*.fasta")]
fasta_finalizer = FastaFinalizer(prime_fasta, final_fasta,
                                 args.p.split())
pool = Pool(processes=args.cpu_number)
pool.map(fasta_finalizer.final_fasta, fasta)
pool.close()
shutil.rmtree(prime_fasta)


# # report generation
# # create database from all tarean_reports
# db_constructor = rtc()
# clusters_table = pd.DataFrame()
# for path, prefix in work_dirs.items():
#     repex_path = path + "/cluster_report.html"
#     table = db_constructor.process_cluster_data(prefix,
#                                                 repex_path,
#                                                 str(args.out) + "/report_src/graph_layouts/" + prefix)
#     clusters_table = clusters_table.append(table, ignore_index=True)
# recomp_results_table = db_constructor.recomp_results_database_construct(str(fasta_path_id),
#                                                                         str(fasta_path_not_id),
#                                                                         str(unique_prob_scl),
#                                                                         str(ref_superclusters),
#                                                                         args.known)
# report_table = db_constructor.recomp_report_table_generation(recomp_results_table, clusters_table)
# report_table.to_csv(str(Path(report_src).parent.parent) + "/superclusters_table.csv",
#                     index=False)
# report_generator = hrg(str(Path(report_src).parent.parent) + "/superclusters_table.csv",
#                        str(out_path) + "/report.html")
# report_generator.generate_report()
# shutil.copyfile("report/style1.css", str(out_path) + "/style1.css")
# logging.info("DONE")
