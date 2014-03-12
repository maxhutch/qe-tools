#!/shared/apps/rhel-6.2/tools/python-3.3.3/bin/python3

from sys import argv, exit
from os import system, chdir, mkdir, getcwd
from os import environ
import os.path 
from shutil import copytree, copy2, rmtree
import time
from optparse import OptionParser
import random
import string 


single_core_rankfile = """
rank 0={0} slot=0:0
"""

dual_core_rankfile = """
rank 0={0} slot=0:0

rank 1={0} slot=1:0
"""


single_socket_rankfile = """
rank 0={0} slot=0:0
rank 1={0} slot=0:1
rank 2={0} slot=0:2
rank 3={0} slot=0:3
rank 4={0} slot=0:4
rank 5={0} slot=0:5
rank 6={0} slot=0:6
rank 7={0} slot=0:7
"""

dual_socket_rankfile = """
rank 0={0} slot=0:0
rank 1={0} slot=0:1
rank 2={0} slot=0:2
rank 3={0} slot=0:3
rank 4={0} slot=0:4
rank 5={0} slot=0:5
rank 6={0} slot=0:6
rank 7={0} slot=0:7

rank 8={0} slot=1:0
rank 9={0} slot=1:1
rank 10={0} slot=1:2
rank 11={0} slot=1:3
rank 12={0} slot=1:4
rank 13={0} slot=1:5
rank 14={0} slot=1:6
rank 15={0} slot=1:7
"""

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

with open(environ["PBS_NODEFILE"], 'r') as f:
  host = (f.readline().split())[0]

with open('rankfile', 'w') as f:
  if opts.np == 16:
    f.write(dual_socket_rankfile.format(host))
  elif opts.np == 8:
    f.write(single_socket_rankfile.format(host))
  elif opts.np == 2:
    f.write(dual_core_rankfile.format(host))
  elif opts.np == 1:
    f.write(single_core_rankfile.format(host))
  else:
    print("Don't know how to run " + str(opts.np) + " ranks") 
    f.write(dual_socket_rankfile)

system("mpirun -np "+ str(opts.np) + " --rankfile rankfile run_with_mps.sh " + opts.bindir + "/vasp 2> "+opts.prefix+".err | tee " + opts.prefix + ".out")

#system("mpirun -np "+ str(opts.np) + " " + opts.bindir + "/vasp 2> "+opts.prefix+".err | tee " + opts.prefix + ".out")
#system("mpirun -np "+ str(opts.np) + " --rankfile rankfile " + opts.bindir + "/vasp 2> "+opts.prefix+".err | tee " + opts.prefix + ".out")


