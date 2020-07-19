import os
import subprocess
from pathlib import Path

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


    def check_blast(self, path):
        if "blast/bin" not in os.environ["PATH"]:
            raise Exception(f"\nblast not detecetd in PATH variable\n")
        try:
            subprocess.check_call(["blastn", "-h"], stdout=open(os.devnull, "w"))
        except subprocess.CalledProcessError:
            raise Exception("\nCan't run blastn. Check PATH variable.\n")
        else:
            if float(".".join(subprocess.check_output(["blastn", "-version"]).split()[1].decode("utf-8").split(".")[:2])) < 2.1:
                raise Exception(f"\nblast 2.10.0+ or higher is require\n")
