"""
All configuration for REcomp
"""

PIPELINE_VERSION = "2.1.2-2"

# default run parameters
EVALUE = 1e-05
CHUNK_SIZE = 10000

# consensuses files
CONSENSUS_FILES = {"RANK1": "TAREAN_consensus_rank_1.fasta",
                   "RANK2": "TAREAN_consensus_rank_2.fasta",
                   "RANK3": "TAREAN_consensus_rank_3.fasta",
                   "RANK4": "TAREAN_consensus_rank_4.fasta"}


# parameters for test run
INPUT_DIRS = {"sample1": "test_data/S1",
              "sample2": "test_data/S2"}
PREFIXES = {"sample1": "S1",
            "sample2": "S2"}
REFERENCES = "test_data/test_references.fasta"
OUTPUT_DIR_IO = "test_data/test_output_io"
OUTPUT_DIR_NIO = "test_data/test_output_nio"
CPU_COUNT = "4"
