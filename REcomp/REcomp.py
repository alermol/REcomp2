#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import itertools
import logging
import os
import shutil
import time
from multiprocessing import Pool, cpu_count
from pathlib import Path

import numpy as np
import pandas as pd
import tqdm
from Bio import SeqIO
from sklearn.cluster import AgglomerativeClustering, KMeans

import config
from common.align_fasta import FastaAligner
from common.check_input import CheckInput
from common.prepare_fasta import PrepFasta as pf
from common.prime_fasta import PrimeFastaWriter
from common.prime_fasta_processing import FastaFinalizer
from common.quick_union import QuickUnion
from report.report_generator import HtmlReportGenerator
from report.summary_table import ReportTableConstructor

parser = argparse.ArgumentParser(description=(
    "REcomp2 - pipeline for comparative analysis of potentially unlimited"
    " number of RepeatExplorer results"
),
    epilog="Please report about all bugs")
parser.add_argument("-v", "--version",
                    help="show version",
                    action="version",
                    version=f"REcomp {config.PIPELINE_VERSION}")
parser.add_argument("i", help="path(s) to RE results (top level)",
                    type=str,
                    metavar="path")
parser.add_argument("p", help="prefix(es) for each paths",
                    type=str,
                    metavar="prefix")
parser.add_argument("out", help="path to output directory")
parser.add_argument("-r",
                    "--references",
                    help="path to fasta with references repeats",
                    metavar="REF")
parser.add_argument("-l", help="save logfile in output directory",
                    action="store_true", dest="log")
parser.add_argument("-c", help="number of CPU to use",
                    type=int, default=cpu_count(),
                    dest="cpu_number", metavar="CPU")
parser.add_argument("-io", "--include-other",
                    help=(
                        "include `other` contigs and clusters "
                        "in analysis (default: False)"
                    ),
                    action="store_true", dest="include_other")
parser.add_argument("-ir", "--include-ribosomal",
                    action="store_true",
                    help=(
                        "include rDNA clusters (rank 4) in analysis "
                        "(default: False)"
                    ),
                    dest="include_ribosomal")
parser.add_argument("--evalue",
                    help=(
                        "evalue threshold for alignments for supercluster "
                        "assembly (default: 1e-05)"
                    ),
                    default=config.EVALUE,
                    type=float)
parser.add_argument("--low-memory",
                    help=("use small amount of RAM for 'all to all' "
                          "blast by using small chunk size (1000) but it "
                          "can take much time (default chunk size: 10000)"
                          ),
                    action="store_true",
                    dest="low_memory")
args = parser.parse_args()

# catch assertions
assert len(args.p.split()) == len(
    set(args.p.split())), ("Prefixes are not unique")
assert len(args.i.split()) == len(
    set(args.i.split())), ("Paths are not unique")
assert 0 < args.cpu_number < cpu_count(), ("CPU count is not valid")
assert args.evalue >= 0.0, ("Wrong E-value thershold")

out_path = Path(args.out)
out_path.mkdir(parents=True, exist_ok=True)

# logging
if args.log:
    logging.basicConfig(level=logging.DEBUG,
                        filename=Path(args.out).joinpath("REcomp.log"),
                        format=("\n%(asctime)s - %(funcName)s - "
                                "%(levelname)s -\n%(message)s\n"),
                        filemode="w")
    console = logging.StreamHandler()
    console.setLevel(logging.DEBUG)
    formatter = logging.Formatter(
        "\n%(asctime)s - %(funcName)s - %(levelname)s -\n%(message)s\n")
    console.setFormatter(formatter)
    logging.getLogger("").addHandler(console)
else:
    logging.basicConfig(level=logging.INFO,
                        format=("\n%(asctime)s - %(funcName)s - "
                                "%(levelname)s -\n%(message)s\n"),)


logging.info(
    (
        f"------------------------------------------------------------------\n"
        f"PIPELINE VERSION               :         {config.PIPELINE_VERSION}\n"
        f"                                                                  \n"
        f"AUTHOR                         :    Aleksey Ermolaev              \n"
        f"------------------------------------------------------------------\n"
    )
)
logging.info(args)

# check input
check_input = CheckInput()
check_input.check_blast(os.environ["PATH"])
if args.references is not None:
    check_input.check_references(args.references)
work_dirs = {path: prefix for path, prefix in zip(
    args.i.split(), args.p.split())}
check_table = check_input.print_check_table(work_dirs)
logging.info((f"Pipeline will be in progress in 30 seconds\n"
              f"Check matching of paths to RE results and their "
              f"prefixes\n{check_table}"))
time.sleep(30)


# create folder structure
logging.info("creating directory structure")
fasta_path = Path(args.out).joinpath("fasta")
fasta_path.mkdir(parents=True, exist_ok=True)

prime_fasta = Path(args.out).joinpath("results", "prime_fasta")
prime_fasta.mkdir(parents=True, exist_ok=True)

final_fasta = Path(args.out).joinpath("results", "final_fasta")
final_fasta.mkdir(parents=True, exist_ok=True)


# prepare fasta files with ranks and "others"
logging.info("creating fasta containing all sequences for analysis")
for path, prefix in work_dirs.items():
    fasta_prep = pf(path, args.references, prefix)
    fasta_prep.create_united_fasta(fasta_path,
                                   include_other=args.include_other,
                                   include_ribosomal=args.include_ribosomal)
if args.references:
    with open(Path(fasta_path).joinpath("fasta.fasta"), "a") as fasta:
        for record in SeqIO.parse(args.references, "fasta"):
            SeqIO.write(record, fasta, "fasta")


# chunk fasta for parallel
records_number = 0
record_iter = SeqIO.parse(
    open(Path(fasta_path).joinpath("fasta.fasta")), "fasta")
chunk_size = config.CHUNK_SIZE
if args.low_memory:
    chunk_size = config.CHUNK_SIZE / 10
logging.info(f"chunk size: {int(chunk_size)}")
time.sleep(5)
for i, batch in enumerate(fasta_prep.batch_iterator(record_iter, chunk_size)):
    records_number += len(batch)
    filename = Path(fasta_path).joinpath(f"fasta{i}.fasta")
    with open(filename, "w") as handle:
        count = SeqIO.write(batch, handle, "fasta")
    logging.info(f"saving chunk {'/'.join(filename.parts[-3:])}")

# prepare connectivity table
logging.info("running all to all blast")
fasta_aligner = FastaAligner(args.evalue)
files = [path for path in fasta_path.rglob("*.fasta")
         if any(map(str.isdigit, Path(path).stem))]
pairs = [list(i) for i in itertools.combinations_with_replacement(files, 2)]
print(f"Running in {args.cpu_number} cpu(s) in parallel")
time.sleep(5)
pool = Pool(processes=args.cpu_number)
result = tqdm.tqdm(pool.imap_unordered(fasta_aligner.align_fasta, pairs),
                   total=len(pairs))
blast_table = pd.concat(result)
pool.close()
logging.info("all to all blast finished")
blast_table = blast_table[blast_table["qseqid"] != blast_table["sseqid"]]
logging.info("removing of junk alignments")
if args.include_other:
    kmeans = KMeans(n_clusters=2).fit(blast_table[["qcovs"]].to_numpy())
else:
    kmeans = AgglomerativeClustering(linkage="single").fit(
        blast_table[["qcovs"]].to_numpy())
bt_kmeans = np.concatenate((blast_table.to_numpy(),
                            kmeans.labels_.reshape(-1, 1)), axis=1)
blast_table = pd.DataFrame(data=bt_kmeans[0:, 0:])
blast_table = blast_table.rename(columns={0: "qseqid",
                                          1: "sseqid",
                                          2: "pident",
                                          3: "qcovs",
                                          4: "cluster"})
if (min(list(itertools.chain(*blast_table[blast_table["cluster"] ==
                                          0][["qcovs"]].values.tolist()))) >
    min(list(itertools.chain(*blast_table[blast_table["cluster"] ==
                                          1][["qcovs"]].values.tolist())))):
    ok_cluster = 0
else:
    ok_cluster = 1
ok_qcovs = min(list(itertools.chain(*blast_table[blast_table["cluster"] ==
                                                 ok_cluster][["qcovs"]].
                                    values.tolist())))
ok_pident = min(list(itertools.chain(*blast_table[blast_table["cluster"] ==
                                                  ok_cluster][["pident"]].
                                     values.tolist())))
blast_table = blast_table[blast_table["cluster"] == ok_cluster]
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
quick_union = QuickUnion(len(map_dict))
for pair in blast_table.itertuples(index=False, name=None):
    quick_union.union(pair[0], pair[1])
uf_repr = [int(i) for i in str(quick_union).split()]
cc_num = set(uf_repr)
logging.info(f"{len(cc_num)} superclusters detected")
time.sleep(5)


# create prime fasta with superclusters
if args.include_other:
    logging.info("prepare primary fasta")
else:
    logging.info("prepare final fasta")
rev_map_dict = {v: k for k, v in map_dict.items()}  # invert map dict
prime_fw = PrimeFastaWriter(fasta_path.joinpath("fasta.fasta"),
                            prime_fasta,
                            final_fasta,
                            rev_map_dict,
                            uf_repr)
print(f"Running in {args.cpu_number} cpu(s) in parallel")
time.sleep(5)
pool = Pool(processes=args.cpu_number)
for _ in tqdm.tqdm(pool.imap_unordered(prime_fw.write_fasta, cc_num),
                   total=len(cc_num)):
    pass
pool.close()


# process prime fasta into final
if args.include_other:
    logging.info(
        "cleaning the primary fasta files from excessive 'other' clusters")
    fasta = [path for path in prime_fasta.rglob("*.fasta")]
    fasta_finalizer = FastaFinalizer(prime_fasta, final_fasta,
                                     args.p.split())
    pool = Pool(processes=args.cpu_number)
    pool.map(fasta_finalizer.final_fasta, fasta)
    pool.close()
shutil.rmtree(prime_fasta)


# report generation
# create database from all tarean_reports
db_constructor = ReportTableConstructor()
clusters_table = pd.DataFrame()
for path, prefix in work_dirs.items():
    repex_path = Path(path).joinpath("cluster_report.html")
    recomp_result_path = out_path.joinpath("report",
                                           "graph_layouts",
                                           prefix)
    table = db_constructor.process_cluster_data(prefix,
                                                repex_path,
                                                recomp_result_path)
    clusters_table = clusters_table.append(table, ignore_index=True)
recomp_results_table = (
    db_constructor.recomp_results_database_construct(final_fasta,
                                                     args.references)
)
report_table = (
    db_constructor.recomp_report_table_generation(recomp_results_table,
                                                  clusters_table)
)
report_table.to_csv(out_path.joinpath("report", "superclusters_table.csv"),
                    index=False)
report_generator = (
    HtmlReportGenerator(out_path.joinpath("report", "superclusters_table.csv"),
                        out_path.joinpath("report.html"),
                        args,
                        ok_qcovs,
                        ok_pident)
)
report_generator.generate_report()
try:
    shutil.copyfile(Path.cwd().joinpath("REcomp2",
                                        "REcomp",
                                        "report",
                                        "style1.css"),
                    out_path.joinpath("style1.css"))
except FileNotFoundError:
    shutil.copyfile(Path.cwd().joinpath("report", "style1.css"),
                    out_path.joinpath("style1.css"))
logging.info("DONE")
