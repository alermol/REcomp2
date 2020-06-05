import itertools
import logging
import re
import shutil
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from bs4 import BeautifulSoup as bs
from sqlalchemy import create_engine


class ReportTableConstructor:
    """
    Class contain methods for building of report table based on results of REcomp
    """
    # def __init__ (self):
    #     self.LOGGER = logging.getLogger(__name__)
    #     self.LOGGER.setLevel(logging.DEBUG)

    def parse_tarean_report(self, path):
        with open(path, "r") as tarean_report:
            soup = bs(tarean_report, features="lxml")
            cluster_list = str(soup.find("script", {"type": "application/json"}))
            cluster_list = re.findall("\[\[.*\]\]", cluster_list)[0][1:-1].split("],[")
            cluster_list = [col.strip("[").strip("]") for col in cluster_list]
            cluster_list = [col.split(",") for col in cluster_list]
            return cluster_list


    def __return_graph_layouts_paths(self,
                                     value,
                                     repex_results_path,
                                     recomp_results_path):
        """
        Function does copy of all graph_layout.png from all folders with results
        was analyzed and returns path to copied files
        """
        gl_repex_path = re.search(r"seqclust.*png",
                                  value.split(">")[0]).group(0)
        copy_from_path = str(Path(repex_results_path).parent) + "/" + gl_repex_path
        copy_to_path = str(recomp_results_path) + re.search(r"/dir_CL[0-9]+.*",
                                                            gl_repex_path).group(0)
        Path(copy_to_path).parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(copy_from_path, copy_to_path)
        return copy_to_path


    def process_cluster_data(self,
                             prefix,
                             repex_results_path,
                             recomp_results_path):
        """
        Function does parse cluster_report.html from each folder with RE results
        and return table with results
        """
        cluster_list = self.parse_tarean_report(repex_results_path)
        cluster_num = [re.search(r"\d+", number).group(0)
                       for number in cluster_list[0]]
        cluster_proportion = [float(value) for value in cluster_list[3]]
        cluster_number_of_reads = [int(value) for value in cluster_list[5]]
        cluster_graph_layout = [self.__return_graph_layouts_paths(value,
                                                                  repex_results_path,
                                                                  recomp_results_path)
                                for value in cluster_list[6]]
        cluster_graph_layout = [re.search(r"report.*", path).group(0)
                                for path in cluster_graph_layout]
        cluster_annotation = [re.search(r"[\w\s\(\)]+", value).group(0)
                              for value in cluster_list[10]]
        prefix = [prefix for _ in range(len(cluster_num))]
        cluster_dataframe = pd.DataFrame({"Prefix": prefix,
                                          "Cluster": cluster_num,
                                          "Proportion": cluster_proportion,
                                          "Number_of_reads": cluster_number_of_reads,
                                          "Graph_layout": cluster_graph_layout,
                                          "TAREAN_annotation": cluster_annotation})
        cluster_dataframe = cluster_dataframe.astype({"Prefix": str,
                                                      "Cluster": int,
                                                      "Proportion": float,
                                                      "Number_of_reads": int,
                                                      "Graph_layout": str,
                                                      "TAREAN_annotation": str})
        return cluster_dataframe

    def recomp_results_database_construct(self,
                                          path_to_identified,
                                          path_to_not_identified,
                                          path_to_probable,
                                          path_to_refscl,
                                          path_to_references):
        """
        Function does compile table with data about each superclusters from
        results of REcomp
        """
        references_id = list(SeqIO.to_dict(SeqIO.parse(path_to_references,
                                                       "fasta")).keys())
        paths = [path_to_identified, path_to_not_identified,
                 path_to_probable, path_to_refscl]
        recomp_results_table = pd.DataFrame()
        for path in paths:
            for file in Path(path).glob("*.fasta"):
                scl_records = list(SeqIO.to_dict(SeqIO.parse(file, "fasta")).values())
                records_id = [scl_id.id for scl_id in scl_records]
                records_seq = [scl_seq.seq for scl_seq in scl_records]
                records_source = ["reference" if record in references_id else
                                  "REother_contig" if "Contig" in record and
                                  "refscl" not in record else
                                  "REreference_scl" if "refscl" in record else
                                  "REconsensus" for record in records_id]
                scl_name = [file.stem for _ in range(len(records_id))]
                uniqueness = ["Truly unique" if len(scl_records) == 1 else ""
                              for _ in range(len(scl_records))]
                path_to_fasta = [re.search(r"[a-z\_]+\/[a-z\_]+\/[a-z0-9\_]+\.[a-z\_]+$", str(file)).group(0)
                                 for _ in range(len(records_id))]
                if path == path_to_identified:
                    cluster_type = ["identified" for _ in range(len(scl_name))]
                if path == path_to_not_identified:
                    cluster_type = ["not_identified" for _ in range(len(scl_name))]
                if path == path_to_probable:
                    cluster_type = ["probable_unique" for _ in range(len(scl_name))]
                if path == path_to_refscl:
                    cluster_type = ["reference_scl" for _ in range(len(scl_name))]
                dataframe = pd.DataFrame({"RecordID": records_id,
                                          "RecordSeq": records_seq,
                                          "ClusterName": scl_name,
                                          "ClusterType": cluster_type,
                                          "RecordSource": records_source,
                                          "Path_to_fasta": path_to_fasta,
                                          "Uniqueness": uniqueness})
                dataframe = dataframe.astype({"RecordID": str,
                                              "RecordSeq": str,
                                              "ClusterName": str,
                                              "ClusterType": str,
                                              "RecordSource": str,
                                              "Path_to_fasta": str,
                                              "Uniqueness": str})
                recomp_results_table = recomp_results_table.append(dataframe,
                                                                   ignore_index=True)
        return recomp_results_table


    def recomp_report_table_generation(self, recomp_table, repex_table):
        """
        Function creates report table with information about all records
        from REcomp results - superclusters_table.csv
        """
        report_table = pd.DataFrame()
        engine = create_engine("sqlite://", echo=False)
        repex_table.to_sql("repex_table", con=engine, index=False)
        for row in recomp_table.itertuples(index=False):
            if row[4] == "reference":
                table = pd.DataFrame({"Prefix": [""],
                                      "Cluster": [""],
                                      "Proportion": [""],
                                      "Number_of_reads": [""],
                                      "Graph_layout": [""],
                                      "TAREAN_annotation": [""],
                                      "RecordID": [row[0]],
                                      "RecordSeq": [row[1]],
                                      "SuperclusterName": [row[2]],
                                      "SuperclusterType": [row[3]],
                                      "RecordSource": [row[4]],
                                      "Features": ["Reference"],
                                      "Path_to_fasta": [row[5]]})
                report_table = report_table.append(table, ignore_index=True)
            elif row[4] == "REconsensus":
                prefix = re.search(r"^[a-zA-Z0-9]+", row[0]).group(0)
                cluster = re.search(r"[0-9]+", row[0].split("_")[1]).group(0)
                proportion = engine.execute(f"SELECT Proportion \
                    FROM repex_table \
                        WHERE Prefix='{prefix}' \
                            AND Cluster='{cluster}'").fetchall()
                number_of_reads = engine.execute(f"SELECT Number_of_reads \
                    FROM repex_table \
                        WHERE Prefix='{prefix}' \
                            AND Cluster='{cluster}'").fetchall()
                grapth_layout = engine.execute(f"SELECT Graph_layout \
                    FROM repex_table \
                        WHERE Prefix='{prefix}' \
                            AND Cluster='{cluster}'").fetchall()
                tarean_annotation = engine.execute(f"SELECT TAREAN_annotation \
                    FROM repex_table \
                        WHERE Prefix='{prefix}' \
                            AND Cluster='{cluster}'").fetchall()
                table = pd.DataFrame({"Prefix": [prefix],
                                      "Cluster": [cluster],
                                      "Proportion": [proportion[0][0]],
                                      "Number_of_reads": [number_of_reads[0][0]],
                                      "Graph_layout": [grapth_layout[0][0]],
                                      "TAREAN_annotation": [tarean_annotation[0][0]],
                                      "RecordID": [row[0]],
                                      "RecordSeq": [row[1]],
                                      "SuperclusterName": [row[2]],
                                      "SuperclusterType": [row[3]],
                                      "RecordSource": [row[4]],
                                      "Features": [row[6]],
                                      "Path_to_fasta": [row[5]]})
                report_table = report_table.append(table, ignore_index=True)
            elif row[4] == "REreference_scl":
                prefix = re.search(r"^[a-zA-Z0-9]+", row[0]).group(0)
                cluster = re.search(r"[0-9]+",
                                    row[0].split("_")[1].split("Contig")[0]).group(0)
                record_id = "_".join(row[0].split("_")[0:2])
                feature = row[0].split("_")[-2]
                proportion = engine.execute(f"SELECT Proportion \
                    FROM repex_table \
                        WHERE Prefix='{prefix}' \
                            AND Cluster='{cluster}'").fetchall()
                number_of_reads = engine.execute(f"SELECT Number_of_reads \
                    FROM repex_table \
                        WHERE Prefix='{prefix}' \
                            AND Cluster='{cluster}'").fetchall()
                grapth_layout = engine.execute(f"SELECT Graph_layout \
                    FROM repex_table \
                        WHERE Prefix='{prefix}' \
                            AND Cluster='{cluster}'").fetchall()
                tarean_annotation = engine.execute(f"SELECT TAREAN_annotation \
                    FROM repex_table \
                        WHERE Prefix='{prefix}' \
                            AND Cluster='{cluster}'").fetchall()
                table = pd.DataFrame({"Prefix": [prefix],
                                      "Cluster": [cluster],
                                      "Proportion": [proportion[0][0]],
                                      "Number_of_reads": [number_of_reads[0][0]],
                                      "Graph_layout": [grapth_layout[0][0]],
                                      "TAREAN_annotation": [tarean_annotation[0][0]],
                                      "RecordID": [record_id],
                                      "RecordSeq": [row[1]],
                                      "SuperclusterName": [row[2]],
                                      "SuperclusterType": [row[3]],
                                      "RecordSource": [row[4]],
                                      "Features": [feature],
                                      "Path_to_fasta": [row[5]]})
                report_table = report_table.append(table, ignore_index=True)
            else:
                prefix = re.search(r"^[a-zA-Z0-9]+", row[0]).group(0)
                cluster = re.search(r"[0-9]+",
                                    row[0].split("_")[1].split("Contig")[0]).group(0)
                record_id = re.search(r"^[a-zA-Z0-9]+_[a-zA-Z0-9]+_[0-9]+bp",
                                      row[0]).group(0)
                feature = row[0].split("_")[-1]
                proportion = engine.execute(f"SELECT Proportion \
                    FROM repex_table \
                        WHERE Prefix='{prefix}' \
                            AND Cluster='{cluster}'").fetchall()
                number_of_reads = engine.execute(f"SELECT Number_of_reads \
                    FROM repex_table \
                        WHERE Prefix='{prefix}' \
                            AND Cluster='{cluster}'").fetchall()
                grapth_layout = engine.execute(f"SELECT Graph_layout \
                    FROM repex_table \
                        WHERE Prefix='{prefix}' \
                            AND Cluster='{cluster}'").fetchall()
                tarean_annotation = engine.execute(f"SELECT TAREAN_annotation \
                    FROM repex_table \
                        WHERE Prefix='{prefix}' \
                            AND Cluster='{cluster}'").fetchall()
                table = pd.DataFrame({"Prefix": [prefix],
                                      "Cluster": [cluster],
                                      "Proportion": [proportion[0][0]],
                                      "Number_of_reads": [number_of_reads[0][0]],
                                      "Graph_layout": [grapth_layout[0][0]],
                                      "TAREAN_annotation": [tarean_annotation[0][0]],
                                      "RecordID": [record_id],
                                      "RecordSeq": [row[1]],
                                      "SuperclusterName": [row[2]],
                                      "SuperclusterType": [row[3]],
                                      "RecordSource": [row[4]],
                                      "Features": [feature],
                                      "Path_to_fasta": [row[5]]})
                report_table = report_table.append(table, ignore_index=True)
        return report_table
