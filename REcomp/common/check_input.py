from Bio import SeqIO
from prettytable import PrettyTable


class CheckInput:

    def print_check_table(self, pp_dict):
        check_table = PrettyTable()
        check_table.field_names = ["Path", "Prefix"]
        for pair in pp_dict.items():
            check_table.add_row([pair[0], pair[1]])
        return check_table


    def check_references(self, known):
        id_list = [record.id for record in SeqIO.parse(known, "fasta")]
        assert len(id_list) == len(set(id_list)), ("Duplicated records in reference FASTA")
