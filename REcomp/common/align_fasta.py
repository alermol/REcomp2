import pandas as pd
from Bio.Blast.Applications import NcbiblastnCommandline


class FastaAligner:

    def __init__(self,
                 evalue,
                 task,
                 database):
        self.evalue = evalue
        self.task = task
        self.database = database

    def align_fasta(self, fasta):
        """
        Create connectivity table for QU using blast
        """
        cline = NcbiblastnCommandline(query=fasta,
                                      db=self.database,
                                      out="-",
                                      outfmt="6 qseqid sseqid pident qcovs",
                                      evalue=self.evalue,
                                      task=self.task)
        output = cline()[0].strip()
        rows = [line.split() for line in output.splitlines()]
        cols = ["qseqid", "sseqid", "pident", "qcovs"]
        data_types = {"qseqid": str,
                      "sseqid": str,
                      "pident": float,
                      "qcovs": float}
        b_tab = pd.DataFrame(rows, columns=cols).astype(data_types)
        return b_tab
