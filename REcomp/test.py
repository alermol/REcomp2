#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Test run
"""
import os

import config
from common import check_input

# check blast in PATH
check_input = check_input.CheckInput()
check_input.check_blast(os.environ["PATH"])

# test run
input_path = (
    f"{' '.join([config.INPUT_DIRS['sample1'], config.INPUT_DIRS['sample2']])}"
)
prefixes = (
    f"{' '.join([config.PREFIXES['sample1'], config.PREFIXES['sample2']])}"
)
# include other
command = (
    f"./REcomp.py '{input_path}' '{prefixes}' {config.OUTPUT_DIR_IO} "
    f" -r {config.REFERENCES} -l -c {config.CPU_COUNT} -io --low-memory"
)
os.system(command)

# not include other
command = (
    f"./REcomp.py '{input_path}' '{prefixes}' {config.OUTPUT_DIR_NIO} "
    f"-r {config.REFERENCES} -l -c {config.CPU_COUNT}"
)
os.system(command)
