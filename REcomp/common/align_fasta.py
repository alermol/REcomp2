import pandas as pd
from Bio.Blast.Applications import NcbiblastnCommandline


class FastaAligner:

    def __init__(self,
                 evalue,
                 perc_identity,
                 qcov_hsp_perc):
        # self.path_to_fasta = path_to_fasta
        self.evalue = evalue
        self.perc_identity = perc_identity
        self.qcov_hsp_perc = qcov_hsp_perc


    def align_fasta(self, pair):
        cline = NcbiblastnCommandline(query=pair[0],
                                      subject=pair[1],
                                      out="-",
                                      outfmt="6 qseqid sseqid",
                                      evalue=self.evalue,
                                      perc_identity=self.perc_identity,
                                      qcov_hsp_perc=self.qcov_hsp_perc)
        print(cline)
        output = cline()[0].strip()
        rows = [line.split() for line in output.splitlines()]
        cols = ["qseqid", "sseqid"]
        data_types = {"qseqid": str, "sseqid": str}
        b_tab = pd.DataFrame(rows, columns=cols).astype(data_types)
        return b_tab
