import pandas as pd
from Bio.Blast.Applications import NcbiblastnCommandline


class FastaAligner:

    def __init__(self, path_to_fasta):
        self.path_to_fasta = path_to_fasta


    def align_fasta(self, fasta):
        cline = NcbiblastnCommandline(query=fasta,
                                      subject=self.path_to_fasta,
                                      out="-",
                                      outfmt="6 qseqid sseqid",
                                      max_hsps=1)
        print(cline)
        output = cline()[0].strip()
        rows = [line.split() for line in output.splitlines()]
        cols = ["qseqid", "sseqid"]
        data_types = {"qseqid": str, "sseqid": str}
        b_tab = pd.DataFrame(rows, columns=cols).astype(data_types)
        return b_tab
