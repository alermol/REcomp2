#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Test run
"""
import os
import subprocess
import sys
from pathlib import Path

import config
from common import check_input

# check blast in PATH
check_input = check_input.CheckInput()
check_input.check_blast(os.environ["PATH"])

# test run
command = (
    f"./REcomp.py '{' '.join([config.INPUT_DIRS['sample1'], config.INPUT_DIRS['sample2']])}' "
    f"'{' '.join([config.PREFIXES['sample1'], config.PREFIXES['sample2']])}' {config.OUTPUT_DIR} "
    f" -r {config.REFERENCES} -l -c {config.CPU_COUNT}"
)
Path(config.OUTPUT_DIR).mkdir(parents=True, exist_ok=True)
os.system(command)
