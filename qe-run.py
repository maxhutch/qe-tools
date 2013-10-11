#!/usr/bin/python3

from sys import argv, exit
from os import system, chdir, mkdir, getcwd
from os.path import exists
from shutil import copytree, copy2, rmtree
import time
from optparse import OptionParser
import random
import string 
import numpy

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
  order = ['control', 'system', 'electrons', 'ions', 'cell', 'shirley']
  for name in order:
    if name in namelists:
      namelists[name].write(f)

def write_cards(cards, f):
  order = ['ATOMIC_SPECIES', 'ATOMIC_POSITIONS', 'K_POINTS', 'Q_POINTS']
  for name in order:
    if name in cards:
      cards[name].write(f)

run_count = 0

def run_pwscf(line, namelists, cards, testdir = '.', np = 1):
  global run_count
  toks = line.split()
  if len(toks) == 2:
    exe = toks[1]
  else:
    exe = 'pw'

  if testdir != '.':
    cwd = getcwd()
    copy2('pw.x', testdir)
    chdir(testdir)

  with open('tmp.in', 'w') as fout:
    write_namelists(namelists, fout)
    write_cards(cards, fout)
 
  system('mpirun -np ' + str(np) + ' ./' + exe + '.x < tmp.in 2> run' + str(run_count) + '.err |tee run' + str(run_count) + '.out') 
  run_count = run_count + 1

  if testdir != '.':
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
parser.add_option("-e", "--executable", dest="exe", default=None,
                  help="Executable directory to be tested")
parser.add_option("-r", "--reference", dest="exe2", default=None,
                  help="Reference executables for comparison testing")
parser.add_option("-n", "--num_pe", type=int, dest="np", default=1,
                  help="Number of MPI PEs")
parser.add_option("-p", "--num_pools", type=int, dest="npool", default=1,
                  help="Number QE Pools")

(opts, args) = parser.parse_args()

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
    run_pwscf(line, namelists, cards, np = opts.np)
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


'''
def run_test(inputs, exe, testdir, np = 1, ipm = False):
  # setup step
#  print("== SETUP =============================================")
  start = time.time()

  rmtree(testdir, ignore_errors = True)
  copytree(inputs, testdir)
  mkdir(testdir + '/bin')
  copy2(exe + '/pw.x', testdir + '/bin')
  copy2(exe + '/shirley_basis.x', testdir + '/bin')
  copy2(exe + '/shirley_ham.x', testdir + '/bin')
  copy2(exe + '/shirley_qdiagp.x', testdir + '/bin')
  copy2(exe + '/bandstruct.x', testdir + '/bin')
  copy2(exe + '../scripts/pwbands.pl', testdir + '/bin')
  copy2(exe + '../scripts/diff_eigvals.pl', testdir + '/bin')
  chdir(testdir)
  time_setup = time.time() - start
#  print(" SETUP completed in %6.2fs" % (time_setup))

#  refloc = './'
  refloc = './ref/'

  # SCF step
  print("-- SCF ----------------------------------------------")
  start = time.time()
  if np == 1:
    system('./bin/pw.x < ' + inputs + '.scf.in > ' + inputs + '.scf.out') 
  else:
    system('mpirun -np ' + str(np) + ' ./bin/pw.x < ' + inputs + '.scf.in > ' + inputs + '.scf.out') 
  end = time.time()
  if exists('CRASH'):
    print('SCF crashed')
    exit(1)
  time_scf = time.time() - start

  with open('./' + inputs + '.scf.out') as out:
    with open(refloc + inputs + '.scf.out') as ref: 
      test_output("!    total energy", 4, "Total Energy", out, ref, 0.001)

  print(" scf completed in %6.2fs" % (time_scf))

  # NSCF step
  print("-- NSCF ---------------------------------------------")
  start = time.time()
  if np == 1:
    system('./bin/pw.x < ' + inputs + '.nscf.in > ' + inputs + '.nscf.out') 
  else:
    system('mpirun -np ' + str(np) + ' ./bin/pw.x < ' + inputs + '.nscf.in > ' + inputs + '.nscf.out') 
  end = time.time()
  if exists('CRASH'):
    print('NSCF crashed')
    exit(1)
  time_nscf = time.time() - start

  with open('./' + inputs + '.nscf.out') as out:
    e_o = sum_energies(out)
  with open(refloc + inputs + '.nscf.out') as ref: 
    e_r = sum_energies(ref)
  compare_vals('Sum E_nk', e_o, e_r)

  print(" nscf completed in %6.2fs" % (time_nscf))

  # Basis step
  print("-- BASIS --------------------------------------------")
  start = time.time()
  if np == 1:
    system('./bin/shirley_basis.x < ' + inputs + '.basis.in > ' + inputs + '.basis.out') 
  else:
    system('mpirun -np ' + str(np) + ' ./bin/shirley_basis.x < ' + inputs + '.basis.in > ' + inputs + '.basis.out') 
  end = time.time()
  if exists('CRASH'):
    print('BASIS crashed')
    exit(1)
  time_basis = time.time() - start

  with open('./' + inputs + '.basis.out') as out:
    with open(refloc + inputs + '.basis.out') as ref: 
      test_output("  -trace(U=-S^2) =  ", 2, "Trace", out, ref, 0.001)
      test_output("  ecutwfc", 2, "Energy Cutoff", out, ref, 0.001)
      test_output(" truncating to", 2, "Basis Size", out, ref, 0.001)

  print(" basis completed in %6.2fs" % (time_basis))

  # Hamiltonian step
  print("-- HAM ----------------------------------------------")
  start = time.time()
  if np == 1:
    system('./bin/shirley_ham.x < ' + inputs + '.ham.in > ' + inputs + '.ham.out') 
  else:
    system('mpirun -np ' + str(np) + ' ./bin/shirley_ham.x < ' + inputs + '.ham.in > ' + inputs + '.ham.out') 
  end = time.time()
  if exists('CRASH'):
    print('HAMIL crashed')
    exit(1)
  time_hamil = time.time() - start
  print(" hamil completed in %6.2fs" % (time_hamil))

  # QDIAG step
  print("-- QDIAG --------------------------------------------")
  start = time.time()
  if np == 1:
    system('./bin/shirley_qdiagp.x < ' + inputs + '.qdiag.in > ' + inputs + '.qdiag.out') 
  else:
    system('mpirun -np ' + str(np) + ' ./bin/shirley_qdiagp.x < ' + inputs + '.qdiag.in > ' + inputs + '.qdiag.out') 
  end = time.time()
  if exists('CRASH'):
    print('QDIAG crashed')
    exit(1)
  time_qdiag = time.time() - start

  with open('./' + inputs + '.qdiag.out') as out:
    with open(refloc + inputs + '.qdiag.out') as ref: 
      test_output("           1", 1, "1st  Eigenval", out, ref, 0.001)
      test_output("        ", 1, "Last Eigenval", out, ref, 0.001)

  print(" qdiag completed in %6.2fs" % (time_qdiag))

  # NSCF check step
  print("-- NSCF2 --------------------------------------------")
  start = time.time()
  if np == 1:
    system('./bin/pw.x < ' + inputs + '.nscf2.in > ' + inputs + '.nscf2.out') 
  else:
    system('mpirun -np ' + str(np) + ' ./bin/pw.x < ' + inputs + '.nscf2.in > ' + inputs + '.nscf2.out') 
  end = time.time()
  if exists('CRASH'):
    print('NSCF2 crashed')
    exit(1)
  time_nscf2 = time.time() - start

  with open('./' + inputs + '.nscf2.out') as out:
    e_o = sum_energies(out)
  with open(refloc + inputs + '.nscf2.out') as ref: 
    e_r = sum_energies(ref)
  compare_vals('Sum E_nk', e_o, e_r)

  print(" nscf2 completed in %6.2fs" % (time_nscf2))

  # Post processing
  print("-- Post Processing ----------------------------------")
  system('./bin/bandstruct.x ' + inputs + '.qdiag.dump > ' + inputs + '.bands.out') 
  system('./bin/pwbands.pl < ' + inputs + '.nscf2.out  > ' + inputs + '.nscf2.eig') 
  system('./bin/diff_eigvals.pl ' + inputs + '.nscf2.eig 5 ' + inputs + '.qdiag.dump.bandstruct 5 > ' + inputs + '.eig.diff')

  chi_o = diff_eigenvalues(inputs + '.nscf2.eig', inputs + '.qdiag.dump.bandstruct')
  chi_r = diff_eigenvalues(refloc + inputs + '.nscf2.eig', refloc + inputs + '.qdiag.dump.bandstruct')
  compare_vals('RMS Error', chi_o, chi_r, 0.001)

  return 1
'''

"""
# read cmdline args
parser = OptionParser()
parser.add_option("-e", "--executable", dest="exe", default=None,
                  help="Executable directory to be tested")
parser.add_option("-r", "--reference", dest="exe2", default=None,
                  help="Reference executables for comparison testing")
parser.add_option("-t", "--test", dest="test", default=None,
                  help="Directory containing the test")
parser.add_option("-s", "--set", dest="test_set", default=None,
                  help="File specifying set of tests to run") 
parser.add_option("-n", "--num_pe", type=int, dest="np", default=1,
                  help="Number of MPI PEs")
parser.add_option("-i", "--ipm", action="store_true", dest="ipm", default=False,
                  help="Use IPM to analyze runs")

(opts, args) = parser.parse_args()

if (opts.test == None and opts.test_set == None) or opts.exe == None:
  parser.error("Must specify at least an executable and test system")

# Setup test set
if opts.test_set == None:
  test_set = [opts.test]
else:
  test_set = []
  with open(opts.test_set) as f:
    for line in f:
      toks = line.split()
      if len(toks) > 0:
        test_set.append(toks[0])

# copy testcase into new director
print("======================================================")    
for test in test_set:

  # Generate run directories 
  testdir = './internal/' + ''.join(random.choice(string.ascii_uppercase + string.digits) for x in range(4))
  if opts.exe2 != None:
    testdir2 = './internal/' + ''.join(random.choice(string.ascii_uppercase + string.digits) for x in range(4))


  print("Test Name: " + test)
  if opts.exe2 == None:
    print("Run on: " + time.strftime("%Y.%m.%d") + "  In: " + testdir)
  else:
   print("Run on: " + time.strftime("%Y.%m.%d") + "  In: " + testdir + "  and: " + testdir2)


  # run first test
  time1 = run_test(test, opts.exe, testdir, opts.np, opts.ipm)

  # if on a comparison, run reference 
  if opts.exe2 != None:
    time2 = run_test(test, opts.exe2, testdir2, opts.np, opts.ipm) 

  print("======================================================")    

"""

'''  
  # Open files
  output = open(testdir + '/OUTCAR')
  if opts.exe2 == None:
    reference = open(testdir + '/EXPECTED')
  else:
    reference = open(testdir2 + '/OUTCAR')

  ref_date = get_tag(reference, "executed on", 4)
  system_name = get_tag(reference, "SYSTEM", 2)
 
  # Test Energy
  print("Result  %-15s %13s vs %13s" % ("Parameter", "Test", "Expected")) 
  print("------------------------------------------------------")    
  test_output("free  energy   TOTEN", 4, "energy", output, reference, .001)
  test_output("external pressure", 3, "ext. pressure", output, reference, .01)
  test_output("volume of cell", 4, "volume", output, reference, .001)
  test_output("  Total  ", 1, "stress (xx)", output, reference, .001)
  test_output("  Total  ", 2, "stress (yy)", output, reference, .001)
  test_output("  Total  ", 3, "stress (zz)", output, reference, .001)
  test_output("  Total  ", 4, "stress (xy)", output, reference, .001)
  test_output("  Total  ", 5, "stress (yz)", output, reference, .001)
  test_output("  Total  ", 6, "stress (zx)", output, reference, .001)
  print("------------------------------------------------------")
  if not opts.ipm:
    compare_output("LOOP:", 6, "loop time", output, reference)
    print("------------------------------------------------------")
    compare_output("SETDIJ:", 6, "setdij time", output, reference)
    compare_output("SUBRT+:", 6, "subrot time", output, reference)
    compare_output("FORHF :", 7, "forhf time", output, reference)
    compare_output("POTLOK:", 6, "potlok time", output, reference)
    compare_output("TRIAL :", 7, "trial time", output, reference)
    compare_output("Elapsed time", 3, "wall time", output, reference)
  print("------------------------------------------------------")
  if opts.exe2 == None:
    print("Tested in %6.2fs" %(time1))
  else:
    print("Test in %6.2fs;  Reference in %6.2fs" % (time1, time2))
  print("======================================================")    

  output.close()
  reference.close()
'''

