import logging
import re
from pathlib import Path

from Bio import SeqIO


class PrepFasta:
    """
    Class contains methods for extracting main fasta files from RE results
    """
    def __init__(self, path_to_results, known_fasta, prefix):
        self.path_to_results = path_to_results
        self.known_fasta = known_fasta
        self.prefix = prefix
        self.RANK1 = "TAREAN_consensus_rank_1.fasta"
        self.RANK2 = "TAREAN_consensus_rank_2.fasta"
        self.RANK3 = "TAREAN_consensus_rank_3.fasta"
        self.RANK4 = "TAREAN_consensus_rank_4.fasta"
        # self.OUT_NAME = prefix + "_ranks123.fasta"
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
            rank_fasta_path = self.path_to_results + rank_file
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
        dirs = Path(self.path_to_results + "/seqclust/clustering/clusters")
        other_dirs = []
        for d in dirs.iterdir():
            if ".." in d.name or d.name in rank_dirs:
                continue
            else:
                other_dirs.append(d.name)
        return other_dirs


    def create_united_fasta(self, path_to_output,
                            exclude_other=True, include_ribosomal=False):
        if not include_ribosomal:
            ranks_files = [self.RANK1, self.RANK2, self.RANK3, self.RANK4]
        else:
            ranks_files = [self.RANK1, self.RANK2, self.RANK3]
        with open(Path(path_to_output).joinpath("fasta.fasta"), "a") as out:
            # ranks
            for file in ranks_files:
                path = self.path_to_results + "/" + file
                with open(path, "r") as rank_fasta:
                    for consensus in SeqIO.parse(rank_fasta, "fasta"):
                        consensus.id = self.prefix + "_" + consensus.id
                        consensus.description = ""
                        SeqIO.write(consensus, out, "fasta")
            # other
            if not exclude_other:
                other_dirs = self.__get_other_dir()
                path = self.path_to_results + "/seqclust/clustering/clusters"
                for file in other_dirs:
                    path_to_cont = path + "/" + file + "/" + "contigs.fasta"
                    with open(path_to_cont, "r") as cont_fasta:
                        for contig in SeqIO.parse(cont_fasta, "fasta"):
                            contig.id = self.prefix + "_" + contig.id
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



    # def get_others_fasta(self, path_to_others):
    #     """
    #     Function does create fasta file(s) of "others" contigs in directory
    #     required by user
    #     """
    #     # logging.info(f"copying of fasta files containing 'other' contigs from {self.path_to_results}")
    #     other_dirs = self.__get_other_dir()
    #     other_fasta_name = path_to_others + "/" + self.prefix + "_other_contigs.fasta"
    #     path = self.path_to_results + "/seqclust/clustering/clusters"
    #     with open(other_fasta_name, "w") as out_fasta:
    #         for od in other_dirs:
    #             path_to_cont = path + "/" + od + "/" + "contigs.fasta"
    #             with open(path_to_cont, "r") as cont_fasta:
    #                 for contig in SeqIO.parse(cont_fasta, "fasta"):
    #                     contig.id = self.prefix + "_" + contig.id
    #                     contig.description = ""
    #                     SeqIO.write(contig, out_fasta, "fasta")


    # def get_ranks_fasta(self, path_to_ranks):
    #     """
    #     Function does create fasta file(s) of ranks consensuses in directory
    #     required by user
    #     """
    #     # logging.info(f"copying of fasta files containing rank's consensuses from {self.path_to_results}")
    #     ranks_files = [self.RANK1, self.RANK2, self.RANK3]
    #     ranks_fasta_path = path_to_ranks + "/" + self.OUT_NAME
    #     for file in ranks_files:
    #         path = self.path_to_results + "/" + file
    #         with open(path, "r") as inp_fasta:
    #             with open(ranks_fasta_path, "a") as out_fasta:
    #                 for consensus in SeqIO.parse(inp_fasta, "fasta"):
    #                     consensus.id = self.prefix + "_" + consensus.id
    #                     consensus.description = ""
    #                     SeqIO.write(consensus, out_fasta, "fasta")
