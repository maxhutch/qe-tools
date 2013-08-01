#!/usr/bin/python3

from sys import argv, exit
from os import system, chdir, mkdir
from os.path import exists
from shutil import copytree, copy2, rmtree
import time
from optparse import OptionParser
import random
import string 
import numpy

# =========================================================================
# Main
#

ref = argv[1]
test = argv[2]
lbnd = int(argv[3])
nbnd = int(argv[4])

if len(argv) > 5:
  efr = float(argv[5])
  eft = float(argv[6])
else:
  efr = 0.
  eft = 0.

if len(argv) > 7:
  pos = int(argv[7])
else:
  pos = 1

ref_dat = numpy.loadtxt(ref)
test_dat = numpy.loadtxt(test)
diff = ref_dat[0:160, 4+lbnd-1:4+nbnd] - test_dat[0:160,4+lbnd-1:4+nbnd] - efr + eft
chi2 = (diff * diff).sum(axis=0) / diff.shape[0]
chit = numpy.sqrt(chi2.sum(axis=0) / len(chi2))
maxloc = numpy.argmax(chi2)
print('Chit: ' + str(chit))
print('Max Chi @ ' + str(maxloc) + ' is ' + str(numpy.sqrt(max(chi2))))


with open('gnu.plt', 'w') as fplt:
  fplt.write("set term wxt \n")
  fplt.write("set nokey \n")
  fplt.write("set style data l \n")
  fplt.write("bandcolor = -1 \n")
  fplt.write('set xtics ("G" 0.0, "H" 0.5, "2H" 1.0)\n')
  fplt.write("set arrow 2 from first 0.5, graph 0 to 0.5, graph 1 lt bandcolor lw 1 nohead \n")
  fplt.write("plot './" + ref + "' u " +str(pos)+ ":(($"+str(4+lbnd)+") - " + str(efr) + ") lt bandcolor \\\n")
  for j in range(5+lbnd,5 + nbnd):
    fplt.write(", '' u " + str(pos)+":(($"+str(j) + ") - " + str(efr) + ") lt bandcolor \\\n")
  fplt.write(", './" + test + "' u " +str(pos)+ ":(($"+str(4+lbnd)+") - " + str(eft) + ") w p \\\n")
  for j in range(5+lbnd,5 + nbnd):
    fplt.write(", '' u "+ str(pos)+":(($"+str(j) + ") - " + str(eft) + ") w p \\\n")
  fplt.write(", 0 \n")

  fplt.write("pause -1 \n")

system('gnuplot gnu.plt')



 



