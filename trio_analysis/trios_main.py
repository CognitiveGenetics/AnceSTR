#!/usr/bin/env python3

import sys, os, subprocess

def main(fileinput, output, build):
    path2script = os.path.join(os.path.dirname(__file__), "trio_analysis", "repeats_inheritance_archive.R")
    path2script = os.path.join(os.path.dirname(__file__), "trio_analysis", "test.R")
    #print(os.path.join(os.path.dirname(__file__), "trio_analysis", "repeats_inheritance.R"))
    #arguments = [fileinput, output, build]
    cmd = ['Rscript', path2script] #+ arguments
    print(cmd)
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    p.stdout.read()