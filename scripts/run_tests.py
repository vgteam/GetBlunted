#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 11:52:19 2020

@author: Jordan
"""

import tempfile
import os
import sys
import subprocess

if __name__ == "__main__":
    
    if len(sys.argv) != 3:
        print("usage: ./run_tests.py bluntify_executable data_dir", file = sys.stderr)
        sys.exit()
    
    here = os.path.dirname(os.path.realpath(__file__))
    validate = os.path.join(here, "validate_bluntification.py")
    bluntify = sys.argv[1]
    datadir = sys.argv[2]
    
    tdir = tempfile.mktemp()
    
    for filename in os.listdir(datadir):
        if not filename.endswith(".gfa") and not filename.endswith(".GFA"):
            continue
        
        in_gfa = os.path.join(datadir, filename)
        print("evaluating GFA {}".format(in_gfa), file = sys.stderr)
        
        prefix = os.path.join(tdir, os.path.basename(filename))
        out_gfa = prefix + ".blunt.gfa"
        out_table = prefix + ".table.txt"
        bluntify_cmd = " ".join([bluntify, in_gfa, out_gfa, out_table])
        print(bluntify_cmd)
        validate_cmd = " ".join([validate, in_gfa, out_gfa, out_table])
        print(validate_cmd)
        #subprocess.check_call(cmd)
        
    