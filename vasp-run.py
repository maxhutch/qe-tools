#!/usr/bin/python3

from sys import argv, exit
from os import system, chdir, mkdir, getcwd
import os.path 
from shutil import copytree, copy2, rmtree
import time
from optparse import OptionParser
import random
import string 

# =========================================================================
# Main
# 

# read cmdline args
parser = OptionParser()
parser.add_option("-e", "--executable", dest="bindir", default='.',
                  help="Executable directory to be tested")
parser.add_option("-r", "--reference", dest="exe2", default=None,
                  help="Reference executables for comparison testing")
parser.add_option("-n", "--num_pe", type=int, dest="np", default=1,
                  help="Number of MPI PEs")
parser.add_option("-p", "--prefix", dest="prefix", default="run",
                  help="Output prefix")

(opts, args) = parser.parse_args()

system("mpirun -np "+ str(opts.np) + " " + opts.bindir + "/vasp 2> "+opts.prefix+".err | tee " + opts.prefix + ".out")

