import itertools
import logging
import re

import numpy as np
import pandas as pd
from natsort import natsorted
from sqlalchemy import create_engine
from yattag import Doc, indent


class HtmlReportGenerator:
    def __init__(self,
                 path_to_report_table,
                 path_to_output_html,
                 args):
        self.path_to_report_table = path_to_report_table
        self.path_to_output_html = path_to_output_html
        self.args = args
        self.LOGGER = logging.getLogger(__name__)
        self.LOGGER.setLevel(logging.DEBUG)


    def generate_report(self):
        logging.info("report generation")
        report_table = pd.read_csv(self.path_to_report_table)
        engine = create_engine("sqlite://", echo=False)
        report_table.to_sql("report_table", con=engine, index=False)
        doc, tag, text = Doc().tagtext()
        doc.asis("<!DOCTYPE html>")
        # header
        with tag("html"):
            with tag("head"):
                doc.stag("link", href="style1.css", rel="stylesheet")
            with tag("h1", klass="main-header"):
                text("REcomp report")
        result = indent(doc.getvalue())
        with open(self.path_to_output_html, "w") as output:
            output.write(result)

        # run parameters
        with tag("h2", klass="run-parameters"):
            text("Run parameters and thresholds:")
            with tag("h4", klass="thresholds"):
                text(f"E-value threshold: {self.args.evalue}")
                doc.stag("br")
                text(f"Identity percent threshold: {self.args.identity_percent}")
                doc.stag("br")
                text(f"Query cover threshold: {self.args.query_cover}")
                doc.stag("br")
                text(f"Include ribosomal clusters: {self.args.include_ribosomal}")
                doc.stag("br")
                text(f"Include 'Other' clusters: {self.args.include_other}")
                doc.stag("br")
        result = indent(doc.getvalue())
        with open(self.path_to_output_html, "w") as output:
            output.write(result)


        # Identified superclusters
        with tag("h2", klass="identified-header"):
            text("Identified superclusters")
        all_identified = engine.execute(f"SELECT SuperclusterName \
            FROM report_table \
                WHERE SuperclusterType=='identified'").fetchall()
        all_identified = list(itertools.chain(*all_identified))
        all_identified = natsorted(list(np.unique(all_identified)))
        for scl in all_identified:
            scl_name = re.search(r"[a-z]+", scl).group(0).capitalize()
            scl_num = re.search(r"\d+", scl).group(0)
            scl_name = f"{scl_name} {scl_num}"
            with tag("h3", klass="supercluster-name"):
                text(scl_name)
            # get path to fasta
            path_to_fasta = engine.execute(f"SELECT Path_to_fasta \
                FROM report_table \
                    WHERE SuperclusterType=='identified' \
                        AND SuperclusterName=='{scl}'").fetchone()[0]
            with tag("a", download="Supercluster fasta", href=path_to_fasta):
                text("Supercluster fasta")
            with tag("table", klass="identified-table"):
                # table header
                with tag("tr"):
                    with tag("th"):
                        text("Cluster")
                    with tag("th"):
                        text("ClusterID")
                    with tag("th"):
                        text("Proportion, %")
                    with tag("th"):
                        text("Number of reads")
                    with tag("th"):
                        text("Sequence length, bp")
                    with tag("th"):
                        text("Sequence")
                    with tag("th"):
                        text("Graph layout")
                    with tag("th"):
                        text("TAREAN annotation")
                    with tag("th"):
                        text("Feature")
                # table for consensuses
                consensuses_id = engine.execute(f"SELECT RecordID FROM report_table \
                    WHERE SuperclusterName=='{scl}' \
                        AND SuperclusterType=='identified' \
                            AND RecordSource=='REconsensus'").fetchall()
                consensuses_id = list(itertools.chain(*consensuses_id))
                for consensus in consensuses_id:
                    scl_info = engine.execute(f"SELECT * FROM report_table \
                        WHERE SuperclusterName=='{scl}' \
                            AND SuperclusterType=='identified' \
                                AND RecordID=='{consensus}'").fetchall()[0]
                    with tag("tr"):
                        with tag("td"):
                            text(int(scl_info[1]))
                        with tag("td"):
                            text(scl_info[6])
                        with tag("td"):
                            text(round(scl_info[2], 3))
                        with tag("td"):
                            text(int(scl_info[3]))
                        with tag("td"):
                            text(len(scl_info[7]))
                        sequence = re.sub("(.{80})", "\\1\n", scl_info[7], 0)
                        with tag("td"):
                            with tag("pre"):
                                text(sequence)
                        with tag("td"):
                            with tag("a", href=scl_info[4]):
                                doc.stag("img",
                                        src=scl_info[4],
                                        width="120", border=0)
                        with tag("td"):
                            text(scl_info[5])
                        with tag("td"):
                            if scl_info[11] is None:
                                text("")
                            else:
                                text(scl_info[11])
                # other contigs
                other_contigs_id = engine.execute(f"SELECT RecordID FROM report_table \
                    WHERE SuperclusterName=='{scl}' \
                        AND SuperclusterType=='identified' \
                            AND RecordSource=='REother_contig'").fetchall()
                other_contigs_id = list(itertools.chain(*other_contigs_id))
                for other in other_contigs_id:
                    scl_info = engine.execute(f"SELECT * FROM report_table \
                        WHERE SuperclusterName=='{scl}' \
                            AND SuperclusterType=='identified' \
                                AND RecordID=='{other}'").fetchall()[0]
                    with tag("tr"):
                        with tag("td"):
                            text(int(scl_info[1]))
                        with tag("td"):
                            text(scl_info[6])
                        with tag("td"):
                            text(round(scl_info[2], 3))
                        with tag("td"):
                            text(int(scl_info[3]))
                        with tag("td"):
                            text(len(scl_info[7]))
                        sequence = re.sub("(.{80})", "\\1\n", scl_info[7], 0)
                        with tag("td"):
                            with tag("pre"):
                                text(sequence)
                        with tag("td"):
                            with tag("a", href=scl_info[4]):
                                doc.stag("img",
                                        src=scl_info[4],
                                        width="120", border=0)
                        with tag("td"):
                            text(scl_info[5])
                        with tag("td"):
                            if scl_info[11] is None:
                                text("")
                            else:
                                text(scl_info[11])

                # references
                references_id = engine.execute(f"SELECT RecordID FROM report_table \
                    WHERE SuperclusterName=='{scl}' \
                        AND SuperclusterType=='identified' \
                            AND RecordSource=='reference'").fetchall()
                references_id = list(itertools.chain(*references_id))
                for reference in references_id:
                    scl_info = engine.execute(f"SELECT * FROM report_table \
                        WHERE SuperclusterName=='{scl}' \
                            AND SuperclusterType=='identified' \
                                AND RecordID=='{reference}'").fetchall()[0]
                    with tag("tr"):
                        with tag("td"):
                            text("")
                        with tag("td"):
                            text(scl_info[6])
                        with tag("td"):
                            text("")
                        with tag("td"):
                            text("")
                        with tag("td"):
                            text(len(scl_info[7]))
                        sequence = re.sub("(.{80})", "\\1\n", scl_info[7], 0)
                        with tag("td"):
                            with tag("pre"):
                                text(sequence)
                        with tag("td"):
                            text("")
                        with tag("td"):
                            text("")
                        with tag("td"):
                            if scl_info[11] is None:
                                text("")
                            else:
                                text(scl_info[11])

        # Not identified superclusters
        with tag("h2", klass="not-identified-header"):
            text("Not identified superclusters")
        not_identified = engine.execute(f"SELECT SuperclusterName \
            FROM report_table \
                WHERE SuperclusterType=='not_identified'").fetchall()
        not_identified = list(itertools.chain(*not_identified))
        not_identified = natsorted(list(np.unique(not_identified)))
        for scl in not_identified:
            scl_name = re.search(r"[a-z]+", scl).group(0).capitalize()
            scl_num = re.search(r"\d+", scl).group(0)
            scl_name = f"{scl_name} {scl_num}"
            with tag("h3", klass="supercluster-name"):
                text(scl_name)
            # get path to fasta
            path_to_fasta = engine.execute(f"SELECT Path_to_fasta \
                FROM report_table \
                    WHERE SuperclusterType=='not_identified' \
                        AND SuperclusterName=='{scl}'").fetchone()[0]
            with tag("a", download="Supercluster fasta", href=path_to_fasta):
                text("Supercluster fasta")
            with tag("table", klass="not-identified-table"):
                # table header
                with tag("tr"):
                    with tag("th"):
                        text("Cluster")
                    with tag("th"):
                        text("ClusterID")
                    with tag("th"):
                        text("Proportion, %")
                    with tag("th"):
                        text("Number of reads")
                    with tag("th"):
                        text("Sequence length, bp")
                    with tag("th"):
                        text("Sequence")
                    with tag("th"):
                        text("Graph layout")
                    with tag("th"):
                        text("TAREAN annotation")
                    with tag("th"):
                        text("Feature")
                # table for consensuses
                consensuses_id = engine.execute(f"SELECT RecordID FROM report_table \
                    WHERE SuperclusterName=='{scl}' \
                        AND SuperclusterType=='not_identified' \
                            AND RecordSource=='REconsensus'").fetchall()
                consensuses_id = list(itertools.chain(*consensuses_id))
                for consensus in consensuses_id:
                    scl_info = engine.execute(f"SELECT * FROM report_table \
                        WHERE SuperclusterName=='{scl}' \
                            AND SuperclusterType=='not_identified' \
                                AND RecordID=='{consensus}'").fetchall()[0]
                    with tag("tr"):
                        with tag("td"):
                            text(int(scl_info[1]))
                        with tag("td"):
                            text(scl_info[6])
                        with tag("td"):
                            text(round(scl_info[2], 3))
                        with tag("td"):
                            text(int(scl_info[3]))
                        with tag("td"):
                            text(len(scl_info[7]))
                        sequence = re.sub("(.{80})", "\\1\n", scl_info[7], 0)
                        with tag("td"):
                            with tag("pre"):
                                text(sequence)
                        with tag("td"):
                            with tag("a", href=scl_info[4]):
                                doc.stag("img",
                                        src=scl_info[4],
                                        width="120", border=0)
                        with tag("td"):
                            text(scl_info[5])
                        with tag("td"):
                            if scl_info[11] is None:
                                text("")
                            else:
                                text(scl_info[11])
                # other contigs
                other_contigs_id = engine.execute(f"SELECT RecordID FROM report_table \
                    WHERE SuperclusterName=='{scl}' \
                        AND SuperclusterType=='not_identified' \
                            AND RecordSource=='REother_contig'").fetchall()
                other_contigs_id = list(itertools.chain(*other_contigs_id))
                for other in other_contigs_id:
                    scl_info = engine.execute(f"SELECT * FROM report_table \
                        WHERE SuperclusterName=='{scl}' \
                            AND SuperclusterType=='not_identified' \
                                AND RecordID=='{other}'").fetchall()[0]
                    with tag("tr"):
                        with tag("td"):
                            text(int(scl_info[1]))
                        with tag("td"):
                            text(scl_info[6])
                        with tag("td"):
                            text(round(scl_info[2], 3))
                        with tag("td"):
                            text(int(scl_info[3]))
                        with tag("td"):
                            text(len(scl_info[7]))
                        sequence = re.sub("(.{80})", "\\1\n", scl_info[7], 0)
                        with tag("td"):
                            with tag("pre"):
                                text(sequence)
                        with tag("td"):
                            with tag("a", href=scl_info[4]):
                                doc.stag("img",
                                        src=scl_info[4],
                                        width="120", border=0)
                        with tag("td"):
                            text(scl_info[5])
                        with tag("td"):
                            if scl_info[11] is None:
                                text("")
                            else:
                                text(scl_info[11])

        # Probable unique superclusters
        with tag("h2", klass="probable-unique-header"):
            text("Probable unique superclusters")
        probable_unique = engine.execute(f"SELECT SuperclusterName \
            FROM report_table \
                WHERE SuperclusterType=='probable_unique' \
                    AND Features!='Truly unique'").fetchall()
        probable_unique = list(itertools.chain(*probable_unique))
        probable_unique = natsorted(list(np.unique(probable_unique)))
        for scl in probable_unique:
            scl_name = re.search(r"[a-z]+", scl).group(0).capitalize()
            scl_num = re.search(r"\d+", scl).group(0)
            scl_name = f"{scl_name} {scl_num}"
            with tag("h3", klass="supercluster-name"):
                text(scl_name)
            # get path to fasta
            path_to_fasta = engine.execute(f"SELECT Path_to_fasta \
                FROM report_table \
                    WHERE SuperclusterType=='probable_unique' \
                        AND SuperclusterName=='{scl}'").fetchone()[0]
            with tag("a", download="Supercluster fasta", href=path_to_fasta):
                text("Supercluster fasta")
            with tag("table", klass="probable-unique-table"):
                # table header
                with tag("tr"):
                    with tag("th"):
                        text("Cluster")
                    with tag("th"):
                        text("ClusterID")
                    with tag("th"):
                        text("Proportion, %")
                    with tag("th"):
                        text("Number of reads")
                    with tag("th"):
                        text("Sequence length, bp")
                    with tag("th"):
                        text("Sequence")
                    with tag("th"):
                        text("Graph layout")
                    with tag("th"):
                        text("TAREAN annotation")
                    with tag("th"):
                        text("Feature")
                # table for consensuses
                consensuses_id = engine.execute(f"SELECT RecordID FROM report_table \
                    WHERE SuperclusterName=='{scl}' \
                        AND SuperclusterType=='probable_unique' \
                            AND RecordSource=='REconsensus'").fetchall()
                consensuses_id = list(itertools.chain(*consensuses_id))
                for consensus in consensuses_id:
                    scl_info = engine.execute(f"SELECT * FROM report_table \
                        WHERE SuperclusterName=='{scl}' \
                            AND SuperclusterType=='probable_unique' \
                                AND RecordID=='{consensus}'").fetchall()[0]
                    with tag("tr"):
                        with tag("td"):
                            text(int(scl_info[1]))
                        with tag("td"):
                            text(scl_info[6])
                        with tag("td"):
                            text(round(scl_info[2], 3))
                        with tag("td"):
                            text(int(scl_info[3]))
                        with tag("td"):
                            text(len(scl_info[7]))
                        sequence = re.sub("(.{80})", "\\1\n", scl_info[7], 0)
                        with tag("td"):
                            with tag("pre"):
                                text(sequence)
                        with tag("td"):
                            with tag("a", href=scl_info[4]):
                                doc.stag("img",
                                        src=scl_info[4],
                                        width="120", border=0)
                        with tag("td"):
                            text(scl_info[5])
                        with tag("td"):
                            if scl_info[11] is None:
                                text("")
                            else:
                                text(scl_info[11])
                # other contigs
                other_contigs_id = engine.execute(f"SELECT RecordID FROM report_table \
                    WHERE SuperclusterName=='{scl}' \
                        AND SuperclusterType=='probable_unique' \
                            AND RecordSource=='REother_contig'").fetchall()
                other_contigs_id = list(itertools.chain(*other_contigs_id))
                for other in other_contigs_id:
                    scl_info = engine.execute(f"SELECT * FROM report_table \
                        WHERE SuperclusterName=='{scl}' \
                            AND SuperclusterType=='probable_unique' \
                                AND RecordID=='{other}'").fetchall()[0]
                    with tag("tr"):
                        with tag("td"):
                            text(int(scl_info[1]))
                        with tag("td"):
                            text(scl_info[6])
                        with tag("td"):
                            text(round(scl_info[2], 3))
                        with tag("td"):
                            text(int(scl_info[3]))
                        with tag("td"):
                            text(len(scl_info[7]))
                        sequence = re.sub("(.{80})", "\\1\n", scl_info[7], 0)
                        with tag("td"):
                            with tag("pre"):
                                text(sequence)
                        with tag("td"):
                            with tag("a", href=scl_info[4]):
                                doc.stag("img",
                                        src=scl_info[4],
                                        width="120", border=0)
                        with tag("td"):
                            text(scl_info[5])
                        with tag("td"):
                            if scl_info[11] is None:
                                text("")
                            else:
                                text(scl_info[11])

        # truly unique
        with tag("h2", klass="truly-unique-header"):
            text("Truly unique superclusters")
        truly_unique = engine.execute(f"SELECT SuperclusterName \
            FROM report_table \
                WHERE SuperclusterType=='probable_unique' \
                    AND Features=='Truly unique'").fetchall()
        truly_unique = list(itertools.chain(*truly_unique))
        truly_unique = natsorted(list(np.unique(truly_unique)))
        with tag("table", klass="probable-unique-table"):
            # table header
            with tag("tr"):
                with tag("th"):
                    text("Cluster")
                with tag("th"):
                    text("ClusterID")
                with tag("th"):
                    text("Proportion, %")
                with tag("th"):
                    text("Number of reads")
                with tag("th"):
                    text("Sequence length, bp")
                with tag("th"):
                    text("Sequence")
                with tag("th"):
                    text("Graph layout")
                with tag("th"):
                    text("TAREAN annotation")
                with tag("th"):
                    text("Feature")
                with tag("th"):
                    text("Fasta")
            # table for consensuses
            for scl in truly_unique:
                consensuses_id = engine.execute(f"SELECT RecordID FROM report_table \
                    WHERE SuperclusterName=='{scl}' \
                        AND SuperclusterType=='probable_unique' \
                            AND RecordSource=='REconsensus'").fetchall()
                consensuses_id = list(itertools.chain(*consensuses_id))
                for consensus in consensuses_id:
                    scl_info = engine.execute(f"SELECT * FROM report_table \
                        WHERE SuperclusterName=='{scl}' \
                            AND SuperclusterType=='probable_unique' \
                                AND RecordID=='{consensus}'").fetchall()[0]
                    with tag("tr"):
                        with tag("td"):
                            text(int(scl_info[1]))
                        with tag("td"):
                            text(scl_info[6])
                        with tag("td"):
                            text(round(scl_info[2], 3))
                        with tag("td"):
                            text(int(scl_info[3]))
                        with tag("td"):
                            text(len(scl_info[7]))
                        sequence = re.sub("(.{80})", "\\1\n", scl_info[7], 0)
                        with tag("td"):
                            with tag("pre"):
                                text(sequence)
                        with tag("td"):
                            with tag("a", href=scl_info[4]):
                                doc.stag("img",
                                        src=scl_info[4],
                                        width="120", border=0)
                        with tag("td"):
                            text(scl_info[5])
                        with tag("td"):
                            if scl_info[11] is None:
                                text("")
                            else:
                                text(scl_info[11])
                        with tag("td"):
                            with tag("a", download="Fasta", href=scl_info[12]):
                                text("FASTA")
        result = indent(doc.getvalue())
        with open(self.path_to_output_html, "w") as output:
            output.write(result)
