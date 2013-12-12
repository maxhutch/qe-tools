#!/usr/bin/python3

from sys import argv, exit
from sys import stdout
from os import system, chdir, mkdir, getcwd
from os import walk
from os.path import exists
from shutil import move, copytree, copy2, rmtree
import time
from optparse import OptionParser
import random
import string 
#import numpy
import json

print_width = 13

print_width = 13

# define test output
def compare_vals(vals, name, tol = 0.001, typ = 'display'):
  output = []
  if vals[0] == None:
    return
  output.append(vals[0])
  for val in vals[1:]:
    if typ == 'intrinsic':
      err = val - vals[0]
    elif typ == 'extrinsic':
      if (vals[0] == 0 and val == 0):
        err = 0.
      elif vals[0] == 0:
        err = 1000000000.
      else:
        err = abs((val - vals[0])/vals[0])
    elif typ == 'display':
      err = val
    elif typ == 'ratio':
      if (vals[0] == 0 and val == 0):
        err = 1.
      else:
        err = val / vals[0]
    output.append(err)

  if typ == 'intrinsic' or typ == 'display':
    fmt = ["{:13s}  {:13.7f}  "]+["{:13.7f}  " if True else "\033[1m{:13.7f}\033[0m  " for x in output[1:]]
#    fmt = ["{:13s}  {:13.7f}  "]+["{:13.7f}  " if abs(x) <= tol else "\033[1m{:13.7f}\033[0m  " for x in output[1:]]
  elif typ == 'extrinsic' or typ == 'ratio':
    fmt = ["{:13s}  {:13.7f}  "]+["{:13.6%}  " if True else "\033[1m{:13.6%}\033[0m  " for x in output[1:]]
#    fmt = ["{:13s}  {:13.7f}  "]+["{:13.6%}  " if abs(x) <= tol else "\033[1m{:13.6%}\033[0m  " for x in output[1:]]
  else:
    print("Compare type "+typ+" not recognized")
  print("".join(fmt).format(name,*output))

# define a value tester
def disp_output(runs, pattern, name, typ = 'display', tol = 1):
  vals = []
  for run in runs:
    val = None
    for out_name in pattern.keys():
      if exists(run+"/"+out_name):
        pat = pattern[out_name][0]
        position = pattern[out_name][1]
        with open(run + "/"+out_name) as f:
          f.seek(0)
          for line in f:
            if pat in line:
              toks = line.split()
              if len(toks) <= position:
                continue
              tmp = toks[position]
              if (tmp[-1] == ':' or 
                  tmp[-1] == "," or 
                  tmp[-1] == ")" or 
                  tmp[-1] == 's' ):
                tmp = tmp[0:-1]
              try: 
                val = float(tmp)
              except ValueError:
                val = 0
        break
    vals.append(val)
  compare_vals(vals, name, tol, typ)
  return vals
'''
def test_eigvals(runs, fermi_energies, tol, nb = -1):   
  if not exists('old/bands.dat'):
    return
  rvals = numpy.loadtxt('old/bands.dat')

  output = ["{:13s}  {:>13s}  ".format("eigvals", " n/a ")]
  for i in range(1,len(runs)):
    run = runs[i]
    if not exists(run + '/bands.dat'):
      print("no " + run + '/bands.dat')
      output.append("{:>13s}  ".format(" n/a ")) 
      continue
    ovals = numpy.loadtxt(run+'/bands.dat')

    nk = min(ovals.shape[0], rvals.shape[0]) 
    if nb < 0:
      nb = min(ovals.shape[1], rvals.shape[1]) - 4

    diff = ovals[0:nk, 4:nb+4] - rvals[0:nk, 4:nb+4]
    if fermi_energies[0] != None:
      diff = diff - fermi_energies[i] + fermi_energies[0]
    chi2 = (diff * diff).sum(axis=1) / nb
    chit = numpy.sqrt(chi2.sum() / nk)

    if chit < tol:
      output.append("{:13.3e}  ".format(chit))
    else:
      output.append("\033[1m{:13.3e}\033[0m  ".format(chit))

  print("".join(output))

  return output
'''
def sum_energies(nscf_out):
  lines = nscf_out.readlines()
  tot = 0.
  for i in range(len(lines)):
    if '          k' in lines[i]:
      for j in range(2,102): # max 100 lines
        toks = lines[i+j].split()
        if len(toks) < 2:
          break
        for k in range(len(toks)):
          tot = tot + float(toks[k])
  return tot

def run_test(inputs, exe, testdir, np = 1, nb = -1, run_mask = None):
  # setup step
  print("------------------------------------------------------")    

  # Copy test into internal directory
  rmtree(testdir, ignore_errors = True)
  copytree(inputs, testdir)
  mkdir(testdir + '/bin')
  if exists(exe + '/vasp'):
    copy2(exe + '/vasp', testdir+'/bin')
  if exists(exe + '/pw.x'):
    copy2(exe + '/pw.x', testdir+'/bin')
  cwd = getcwd()
  chdir(testdir)

  # Setup list of runs
  runs = next(walk('.'))[1];  
  runs.sort();
  runs = [x for x in runs if x != 'pseudo' and x != 'bin']
  # Make a copy of the 'ref' run for comparison
  if 'ref' in runs:
    runs.remove('ref')
    copytree('./ref', './old')
    runs.insert(0, 'old')
    runs.insert(1, 'ref')

  # mask the runs
  if run_mask != None:
    do_runs = run_mask.split(":")
    do_runs.append("old")
    runs = [x for x in runs if x in do_runs]
 
  # Loop over the runs, running each and recording time 
  fmt = "{:13s}  "+"{:>13s}  "*len(runs)
  print(fmt.format("runs", *runs))
  print("------------------------------------------------------")    
  print("{:13s}  ".format("time (s)"), end=""); stdout.flush()
  for run in runs:
    if run == 'old':
      with open('./old/time', 'r') as f:
        print("%13.1f  " % float(f.readline()), end=""); stdout.flush()
      continue
    chdir(run)
    start = time.time()  
    if exists('./INCAR'):
      system('vasp-run.py -n ' + str(np) + ' -e ../bin/ > /dev/null')
    else:
      system('qe-run.py ' + inputs + '.'+run+'.in -n ' + str(np) + ' -e ../bin/ > /dev/null')
    time_test = time.time() - start
    system('echo "'+str(time_test)+'" > time')
    print("%13.1f  " % time_test , end=""); stdout.flush()
    chdir('..')
  print("")

  # set tolerances etc.
  loaded_config = {}
  if exists('./config'):
    with open('./config', 'r') as f:
      loaded_config = json.load(f)
  default_config = {
                    "nb":     nb,
                    "etot":   7.34e-5,
                    "efermi": 0.01,
                    "force":  0.01,
                    "stress": 1.,
                    "bands":  0.01
                   } 
  config = dict(list(default_config.items()) + list(loaded_config.items()))
  total_energy = {"OUTCAR": ["free  energy   TOTEN", 4], "run0.out": ["!    total energy", 4]}
  fermi_energy = {"OUTCAR": ["E-fermi", 2], "run0.out": ["the Fermi energy", 4]}
  total_force = {"run0.out":  ["Total force =", 3]}
  pressure =     {"OUTCAR": ["external pressure", 3], "run0.out": ["total   stress", 5]} 
  # Compare some fields
  disp_output(runs, total_energy, "Total Energy", 'extrinsic', config['etot'])
  disp_output(runs, fermi_energy, "Fermi Energy", 'intrinsic', config['efermi'])
  disp_output(runs, total_force, "Total Force", 'extrinsic', config['force'])
  disp_output(runs, pressure, "Pressure", 'intrinsic', config['stress'])

  print("------------------------------------------------------")    

  gvecs = {"OUTCAR": ["total plane-waves  NPLWV", 4], "run0.out": ["Kohn-Sham Wavefunctions", 5]}
  ecut = {"run0.out": ["kinetic-energy cutoff", 3]}

  disp_output(runs, ecut, "Energy Cutoff")
  disp_output(runs, gvecs, "# Gvectors")

  print("------------------------------------------------------")    
  scf_time = {"OUTCAR": ["LOOP+:", 3], "run0.out": ["electrons    :", 4]}
  nscf_time = {"run1.out": ["electrons    :", 4]}
  forces_time = {"run0.out": ["forces       :", 4]}
  stresses_time = {"run0.out": ["stress       :", 4]}
#  disp_output(runs, scftime, "SCF Time", 'ratio')
  disp_output(runs, scf_time, "SCF Time", 'ratio')
  disp_output(runs, nscf_time, "NSCF Time", 'ratio')
  disp_output(runs, forces_time, "Force Time", 'ratio')
  disp_output(runs, stresses_time, "Stress Time", 'ratio')
  chdir(cwd)
  return 1
'''
  ef = disp_output(runs, "the Fermi energy", 4, "Fermi Energy", config['efermi'])
  disp_output(runs, "convergence has been achieved in", 5, "N Iter", 1)
  disp_output(runs, "Total force", 3, "Total Force", config['force'], intrinsic=False)
  disp_output(runs, "temperature", 2, "Temperature", 0.01)

  test_eigvals(runs, ef, config["bands"], config["nb"])
'''

'''

  move('./run0.out', './run0.test.out')
#  move('./run1.out', './run1.test.out')
  move('./run0.err', './run0.test.err')
#  move('./run1.err', './run1.test.err')
  if compare_bands:
    move('./bands.dat', './bands.test.dat')
  print(" Test ran in      %6.2fs" % (time_test))

  if rerun:
    start = time.time()
    system('qe-run.py ' + inputs + '.ref.in -n ' + str(np) + ' > /dev/null')
    time_ref = time.time() - start
    move('./run0.out', './run0.ref.out')
#    move('./run1.out', './run1.ref.out')
    move('./run0.err', './run0.ref.err')
#    move('./run1.err', './run1.ref.err')
    if compare_bands:
      move('./bands.dat', './bands.ref.dat')
    print(" Reference ran in %6.2fs" % (time_ref))
  else:
    move('ref/run0.ref.out', './run0.ref.out')
#    move('ref/run1.ref.out', './run1.ref.out')
    move('ref/run0.ref.err', './run0.ref.err')
#    move('ref/run1.ref.err', './run1.ref.err')
    if compare_bands:
      move('ref/bands.ref.dat', './bands.ref.dat')
'''

# read cmdline args
parser = OptionParser()
parser.add_option("-e", "--executable", dest="exe", default=None,
                  help="Executable directory to be tested")
parser.add_option("-r", "--runs", dest="runs", default=None,
                  help="Only perform these runs (: separated)")
parser.add_option("-t", "--test", dest="test", default=None,
                  help="Directory containing the test")
parser.add_option("-s", "--set", dest="test_set", default=None,
                  help="File specifying set of tests to run") 
parser.add_option("-n", "--num_pe", type=int, dest="np", default=1,
                  help="Number of MPI PEs")
#parser.add_option("-f", "--force", action="store_true", dest="force", default=False,
#                  help="Force the reference run to occur")
#parser.add_option("-i", "--ipm", action="store_true", dest="ipm", default=False,
#                  help="Use IPM to analyze runs")
parser.add_option("-b", "--nband", type=int, dest="nb", default=-1,
                  help="How many bands to compare")


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
      if line[0] == '!':
        continue 
      toks = line.split()
      if len(toks) > 0:
        test_set.append(toks[0])

# copy testcase into new director
print("======================================================")    
for test in test_set:

  # Generate run directories 
  testdir = './internal/' + ''.join(random.choice(string.ascii_uppercase + string.digits) for x in range(4))

  print("Test Name: " + test)
  print("Run on: " + time.strftime("%Y.%m.%d") + "  In: " + testdir)


  # run first test
  time1 = run_test(test, opts.exe, testdir, opts.np, opts.nb, opts.runs)

  print("======================================================")    

