import pandas as pd
from Bio.Blast.Applications import NcbiblastnCommandline


class FastaAligner:

    def __init__(self,
                 evalue):
        self.evalue = evalue

    def align_fasta(self, pair):
        """
        Create connectivity table for QU using blast
        """
        cline = NcbiblastnCommandline(query=pair[0],
                                      subject=pair[1],
                                      out="-",
                                      outfmt="6 qseqid sseqid pident qcovs",
                                      evalue=self.evalue)
        output = cline()[0].strip()
        rows = [line.split() for line in output.splitlines()]
        cols = ["qseqid", "sseqid", "pident", "qcovs"]
        data_types = {"qseqid": str,
                      "sseqid": str,
                      "pident": float,
                      "qcovs": float}
        b_tab = pd.DataFrame(rows, columns=cols).astype(data_types)
        return b_tab
