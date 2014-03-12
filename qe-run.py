#!/shared/apps/rhel-6.2/tools/python-3.3.3/bin/python3

from sys import argv, exit
from os import system, chdir, mkdir, getcwd
import os.path 
from shutil import copytree, copy2, rmtree
import time
from optparse import OptionParser
import random
import string 

class Namelist:
  def __init__(self, name):
    self.name = name
    self.params = {}

  def write(self, f):
    f.write(' &' + self.name + '\n')
    for param in self.params.items():
      f.write('  ' + param[0] + ' = ' + param[1] + ' \n')
    f.write(' / \n')

  def add(self, name, val):
    self.params[name] = val

class Card:
  def __init__(self, name):
    self.name = name
    self.lines = []

  def add(self, line):
    self.lines.append(line)

  def write(self, f):
    for line in self.lines:
      f.write(line)

def parse_line(line, namelists, cards, active_nl, active_card):
  '''
  if active_nl[0] != None:
    print('NL ' + active_nl[0].name + ' is active')
  if active_card[0] != None:
    print('Card ' + active_card[0].name + ' is active')
  if active_nl[0] == None and active_card[0] == None:
    print('Neither NL nor Card active \n')
  '''
  toks = line.split()
  if len(toks) == 0:
    active_nl[0] = None
    active_card[0] = None
    return

  if toks[0][0] == '!' or toks[0][0] == '@':
    return
  if toks[0] == '/':
    active_nl[0] = None
    return
  if toks[0][0] == '&':
    active_card[0] = None
    if toks[0][1:] in namelists:
      active_nl[0] = namelists[toks[0][1:]]
    else:
      namel = Namelist(toks[0][1:])
      namelists[toks[0][1:]] = namel
      active_nl[0] = namel
    return
  if active_nl[0] != None:
    if toks[1] != '=':
      print('Invalid param: ' + line)
    active_nl[0].add(toks[0], toks[2])
    return
  if active_card[0] != None:
    active_card[0].add(line)
    return  
  if toks[0] in cards:
    del cards[toks[0]]
  new_card = Card(toks[0])
  active_card[0] = new_card
  cards[toks[0]] = new_card
  new_card.add(line)

  return

def write_namelists(namelists, f):
  order = ['control', 'system', 'electrons', 'ions', 'cell', 'srb', 'shirley']
  for name in order:
    if name in namelists:
      namelists[name].write(f)

def write_cards(cards, f):
  order = ['ATOMIC_SPECIES', 'ATOMIC_POSITIONS', 'CELL_PARAMETERS', 'K_POINTS', 'Q_POINTS']
  for name in order:
    if name in cards:
      cards[name].write(f)

run_count = 0

def run_pwscf(line, namelists, cards, opts):
  global run_count
  toks = line.split()
  if len(toks) == 2:
    exe = os.path.join(opts.bindir, toks[1]+'.x')
  else:
    exe = os.path.join(opts.bindir, 'pw.x')

  if opts.testdir != '.':
    cwd = getcwd()
    copy2('pw.x', opts.testdir)
    chdir(opts.testdir)

  with open('tmp.in', 'w') as fout:
    write_namelists(namelists, fout)
    write_cards(cards, fout)

  system('mpirun -np ' + str(opts.nproc) + ' ' + exe + ' -npot ' + str(opts.npot) + ' -npool ' + str(opts.npool) + ' < tmp.in 2> '+ opts.prefix + str(run_count) + '.err |tee '+opts.prefix+ str(run_count) + '.out') 
  run_count = run_count + 1

  if opts.testdir != '.':
    chdir(cwd)

  return

def make_bands(line):
  with open('gnu.plt', 'w') as fplt:
    fplt.write("set nokey \n")
    fplt.write("set style data l \n")
    fplt.write("bandcolor = -1 \n")
    fplt.write('set xtics ("G" 0.0, "H" 0.5, "2H" 1.0)\n')
    fplt.write("set arrow 2 from first 0.5, graph 0 to 0.5, graph 1 lt bandcolor lw 1 nohead \n")
    fplt.write("plot './bands.dat' u 2:5 lt bandcolor \\\n")
    fplt.write(", '' u 2:6 lt bandcolor \\\n")
    fplt.write(", '' u 2:7 lt bandcolor \\\n")
    fplt.write(", '' u 2:8 lt bandcolor \\\n")
    fplt.write(", '' u 2:9 lt bandcolor \\\n")
    fplt.write(", '' u 2:10 lt bandcolor \\\n")
    fplt.write(", '' u 2:11 lt bandcolor \\\n")
    fplt.write(", '' u 2:12 lt bandcolor  \n")
    fplt.write("pause mouse  \n")
  system('gnuplot gnu.plt')

def cbands(line):
  toks = line.split()
  with open('gnu.plt', 'w') as fplt:
    fplt.write("set nokey \n")
    fplt.write("set style data l \n")
    fplt.write("bandcolor = -1 \n")
    fplt.write('set xtics ("G" 0.0, "H" 0.5, "2H" 1.0)\n')
    fplt.write("set arrow 2 from first 0.5, graph 0 to 0.5, graph 1 lt bandcolor lw 1 nohead \n")
    fplt.write("plot './" + toks[1] + "' u 2:5 lt bandcolor \\\n")
    fplt.write(", '' u 2:6 lt bandcolor \\\n")
    fplt.write(", '' u 2:7 lt bandcolor \\\n")
    fplt.write(", '' u 2:8 lt bandcolor \\\n")
    fplt.write(", '' u 2:9 lt bandcolor \\\n")
    fplt.write(", '' u 2:10 lt bandcolor \\\n")
    fplt.write(", '' u 2:11 lt bandcolor \\\n")
    fplt.write(", '' u 2:12 lt bandcolor \\\n")
    fplt.write(", './" + toks[2] + "' u 2:5 w p \\\n")
    fplt.write(", '' u 2:6 w p \\\n")
    fplt.write(", '' u 2:7 w p \\\n")
    fplt.write(", '' u 2:8 w p \\\n")
    fplt.write(", '' u 2:9 w p \\\n")
    fplt.write(", '' u 2:10 w p \\\n")
    fplt.write(", '' u 2:11 w p \\\n")
    fplt.write(", '' u 2:12 w p  \n")

    fplt.write("pause mouse  \n")
  system('gnuplot gnu.plt')



# =========================================================================
# Main
# 

# read cmdline args
parser = OptionParser()
parser.add_option("-e", "--executable", dest="bindir", default='.',
                  help="Executable directory to be tested")
parser.add_option("-r", "--reference", dest="exe2", default=None,
                  help="Reference executables for comparison testing")
parser.add_option("--nproc", type=int, dest="nproc", default=1,
                  help="Number of MPI PEs")
parser.add_option("--npool", type=int, dest="npool", default=1,
                  help="Number of pools")
parser.add_option("--npot", type=int, dest="npot", default=-1,
                  help="Number of pots (srb pools)")
parser.add_option("-p", "--prefix", dest="prefix", default="run",
                  help="Output prefix")
parser.add_option("-d", "--dir", dest="testdir", default=".",
                  help="Directory containing test")



(opts, args) = parser.parse_args()

if opts.npot == -1:
  opts.npot = opts.nproc

fin = open(argv[1], 'r')
namelists = {}
cards = {}
active_nl = [None]
active_card = [None]
lines = fin.readlines()
for line in lines:
  parse_line(line, namelists, cards, active_nl, active_card) 
  toks = line.split()
  if len(toks) == 0 or toks[0][0] != '@': 
    continue
  if toks[0] == '@run':
    run_pwscf(line, namelists, cards, opts)
  elif toks[0] == '@bands':
    make_bands(line)
  elif toks[0] == '@cmd':
    cmdline = ''
    for tok in toks[1:]:
      cmdline = cmdline + tok + ' '
    system(cmdline)
  elif toks[0] == '@cbands':
    cbands(line)
fin.close()

