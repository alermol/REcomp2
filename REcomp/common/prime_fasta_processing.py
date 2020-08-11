import pandas as pd
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline


class FastaFinalizer:

    def __init__(self,
                 path_to_prime,
                 path_to_final,
                 prefix_list):
        self.path_to_prime = path_to_prime
        self.path_to_final = path_to_final
        self.prefix_list = prefix_list

    def final_fasta(self, path_to_fasta):
        records_id = [record.id for record in SeqIO.parse(path_to_fasta,
                                                          "fasta")]
        rank_prefixes = [rank_id.split("_")[0] for rank_id in records_id
                         if "_TR_1_x_" in rank_id]
        other_id = [other for other in records_id if "Contig" in other]
        other_prefixes = [oth_id.split("_")[0] for oth_id in other_id]
        require_pref = set(other_prefixes).difference(set(rank_prefixes))
        ffasta_path = self.path_to_final.joinpath(path_to_fasta.name)
        prime_tmp_fasta_path = self.path_to_prime.joinpath(
            f"{path_to_fasta.stem}_tmp.fasta"
        )
        prime_fasta_path = self.path_to_prime.joinpath(path_to_fasta.name)
        with open(ffasta_path, "w") as ffasta:
            for record in SeqIO.parse(path_to_fasta, "fasta"):
                if record.id not in other_id:
                    SeqIO.write(record, ffasta, "fasta")
                else:
                    continue
        if len(require_pref) == 0:
            return
        tmp_id = self.__tmp_other_id(require_pref, other_id)
        with open(prime_tmp_fasta_path, "w") as tmp, \
                open(prime_fasta_path) as src:
            for record in SeqIO.parse(src, "fasta"):
                if record.id in tmp_id:
                    SeqIO.write(record, tmp, "fasta")
        cline = NcbiblastnCommandline(query=ffasta_path,
                                      subject=prime_tmp_fasta_path,
                                      out="-",
                                      outfmt="6 qseqid sseqid evalue bitscore",
                                      task="blastn")
        output = cline()[0].strip()
        rows = [line.split() for line in output.splitlines()]
        cols = ["qseqid", "sseqid", "evalue", "bitscore"]
        data_types = {"qseqid": str,
                      "sseqid": str,
                      "evalue": float,
                      "bitscore": float}
        b_tab = pd.DataFrame(rows, columns=cols).astype(data_types)
        if b_tab.empty:
            return
        b_tab_bit = (
            b_tab.sort_values("bitscore").iloc[[0]]["sseqid"].tolist()[0]
        )
        b_tab_evalue = (
            b_tab.sort_values("evalue").iloc[[0]]["sseqid"].tolist()[0]
        )
        with open(ffasta_path, "a") as ffasta, \
                open(prime_tmp_fasta_path) as tmp:
            for record in SeqIO.parse(tmp, "fasta"):
                if record.id == b_tab_bit:
                    record.id = f"{record.id}_maxBitscore"
                    record.description = ""
                    SeqIO.write(record, ffasta, "fasta")
                elif record.id == b_tab_evalue:
                    record.id = f"{record.id}_minEvalue"
                    record.description = ""
                    SeqIO.write(record, ffasta, "fasta")
                else:
                    continue

    def __tmp_other_id(self, pref_list, other_id):
        other_id_tmp = []
        for pref in pref_list:
            for oth_id in other_id:
                if pref in oth_id:
                    other_id_tmp.append(oth_id)
        return other_id_tmp
