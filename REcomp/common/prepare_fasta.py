import logging
import re
from pathlib import Path

from Bio import SeqIO

import config


class PrepFasta:
    """
    Class contains methods for extracting main fasta files from RE results
    """

    def __init__(self, path_to_results, known_fasta, prefix):
        self.path_to_results = path_to_results
        self.known_fasta = known_fasta
        self.prefix = prefix
        self.RANK1 = config.CONSENSUS_FILES["RANK1"]
        self.RANK2 = config.CONSENSUS_FILES["RANK2"]
        self.RANK3 = config.CONSENSUS_FILES["RANK3"]
        self.RANK4 = config.CONSENSUS_FILES["RANK4"]
        self.LOGGER = logging.getLogger(__name__)
        self.LOGGER.setLevel(logging.DEBUG)

    def __get_ranks_dir(self):
        """
        Function does generate list of directories names of clusters in
        rank fasta (e.g. convert CL10 into dir_CL0010)
        """
        rank_files = [self.RANK1, self.RANK2, self.RANK3, self.RANK4]
        known_cl = []
        for rank_file in rank_files:
            rank_fasta_path = self.path_to_results + "/" + rank_file
            rank_fasta = SeqIO.to_dict(SeqIO.parse(rank_fasta_path, "fasta"))
            for cl in list(rank_fasta.keys()):
                cl_part = re.search(r"CL\d{1,}", cl).group(0)
                cl_num_part = "000" + re.search(r"\d+", cl_part).group(0)
                known_cl.append("dir_" + "CL" + cl_num_part[-4:])
        return known_cl

    def __get_other_dir(self):
        """
        Function does generate list of directories names of "other" clusters
        """
        rank_dirs = self.__get_ranks_dir()
        dirs = Path(self.path_to_results).joinpath("seqclust",
                                                   "clustering",
                                                   "clusters")
        other_dirs = []
        for d in dirs.iterdir():
            if ".." in d.name or d.name in rank_dirs:
                continue
            else:
                other_dirs.append(d.name)
        return other_dirs

    def create_united_fasta(self, path_to_output,
                            include_other=False, include_ribosomal=False):
        if include_ribosomal:
            ranks_files = [self.RANK1, self.RANK2, self.RANK3, self.RANK4]
        else:
            ranks_files = [self.RANK1, self.RANK2, self.RANK3]
        with open(Path(path_to_output).joinpath("fasta.fasta"), "a") as out:
            # ranks
            for file in ranks_files:
                path = Path(self.path_to_results).joinpath(file)
                with open(path, "r") as rank_fasta:
                    for consensus in SeqIO.parse(rank_fasta, "fasta"):
                        consensus.id = f"{self.prefix}_{consensus.id}"
                        consensus.description = ""
                        SeqIO.write(consensus, out, "fasta")
            # other
            if include_other:
                other_dirs = self.__get_other_dir()
                path = Path(self.path_to_results).joinpath("seqclust",
                                                           "clustering",
                                                           "clusters")
                for file in other_dirs:
                    path_to_cont = Path(path).joinpath(file, "contigs.fasta")
                    with open(path_to_cont, "r") as cont_fasta:
                        for contig in SeqIO.parse(cont_fasta, "fasta"):
                            contig.id = f"{self.prefix}_{contig.id}"
                            contig.description = ""
                            SeqIO.write(contig, out, "fasta")

    def batch_iterator(self, iterator, batch_size):
        entery = True
        while entery:
            batch = []
            while len(batch) < batch_size:
                try:
                    entery = next(iterator)
                except StopIteration:
                    entery = None
                if entery is None:
                    break
                batch.append(entery)
            if batch:
                yield batch
