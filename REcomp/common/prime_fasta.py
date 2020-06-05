import re
from pathlib import Path

from Bio import SeqIO


class PrimeFastaWriter:

    def __init__(self,
                 input_fasta,
                 prime_output,
                 final_output,
                 map_dict,
                 uf_list):
        self.input_fasta = input_fasta
        self.prime_output = prime_output
        self.final_output = final_output
        self.map_dict = map_dict
        self.uf_list = uf_list


    def write_fasta(self, scl_num):
        ind_list = [i for i,value in enumerate(self.uf_list)
                    if value == scl_num]
        id_list = [self.map_dict[i] for i in ind_list]
        if (len(id_list) == 1 and "_TR_1_x_" not in id_list[0] or
            len(id_list) == 1 and "Contig" in id_list[0] or
            not any("_TR_1_x_" in i for i in id_list)):
            return
        else:
            if (len(id_list) == 1 and "_TR_1_x_" in id_list[0] or
                not any("Contig" in i for i in id_list)):
                file = Path(self.final_output).joinpath(f"supercluster{scl_num}.fasta")
            else:
                file = Path(self.prime_output).joinpath(f"supercluster{scl_num}.fasta")
            with open(file, "w") as output_fasta:
                for record in SeqIO.parse(open(self.input_fasta), "fasta"):
                    if record.id in id_list:
                        SeqIO.write(record, output_fasta, "fasta")
