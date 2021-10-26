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
import hashlib
import time

from collections import defaultdict

if __name__ == "__main__":
    
    if len(sys.argv) != 3:
        print("usage: ./run_tests.py bluntify_executable data_dir", file=sys.stderr)
        sys.exit()
    
    here = os.path.dirname(os.path.realpath(__file__))
    validate = os.path.join(here, "validate_bluntification.py")
    generate = os.path.join(here, "generate_gfa.py")
    bluntify = sys.argv[1]
    datadir = sys.argv[2]
    
    tdir = tempfile.mktemp()
    os.makedirs(tdir)

    all_pass = True
    failing_tests = defaultdict(list)

    for filename in os.listdir(datadir):
        if not filename.endswith(".gfa") and not filename.endswith(".GFA"):
            continue
        
        in_gfa_path = os.path.join(datadir, filename)
        print("evaluating GFA {}".format(in_gfa_path), file=sys.stderr)
        
        prefix = os.path.join(tdir, os.path.basename(filename))
        out_gfa_path = prefix + ".blunt.gfa"
        out_table_path = prefix + ".table.txt"

        with open(out_gfa_path, 'w') as out_gfa:
            bluntify_cmd = [bluntify, in_gfa_path, "-p", out_table_path]
            print(" ".join(bluntify_cmd))

            try:
                subprocess.run(bluntify_cmd, stdout=out_gfa, check=True)
            except Exception as e:
                print("FAIL")
                print(e)

                failing_tests[" ".join(bluntify_cmd)].append(e)
                all_pass = False


            validate_cmd = [validate, in_gfa_path, out_gfa_path, out_table_path]
            print(" ".join(validate_cmd))

            try:
                subprocess.run(validate_cmd, stdout=out_gfa, check=True)
            except Exception as e:
                print(e)
                failing_tests[" ".join(validate_cmd)].append(e)

                all_pass = False

    if not all_pass:
        for item in failing_tests.items():
            print("FAIL: ", item[0])
            for subitem in item[1]:
                print(subitem)
            print()

        exit("FAIL")

    # num_random_gfas = 50
    # for i in range(num_random_gfas):
    #     gfa_name = os.path.join(tdir, hashlib.sha256(str.encode(str(time.time()))).hexdigest()[:10] + ".gfa")
    #     generate_cmd = generate + " > " + gfa_name
    #     print("evaluating GFA {}".format(gfa_name), file = sys.stderr)
    #
    #     prefix = os.path.join(tdir, os.path.basename(gfa_name))
    #     out_gfa = prefix + ".blunt.gfa"
    #     out_table = prefix + ".table.txt"
    #     bluntify_cmd = " ".join([bluntify, in_gfa, out_gfa, out_table])
    #     print(bluntify_cmd)
    #     validate_cmd = " ".join([validate, in_gfa, out_gfa, out_table])
    #     print(validate_cmd)
