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
import numpy
import json
from math import sqrt

print_width = 13

print_width = 13

# define test output
def compare_vals(vals, name, tol = 0.001, typ = 'display'):
  output = []
  if vals[0] == None:
    return
  output.append(vals[0])
  fmt = ["{:13s}  {:13.7f}  "]
  for val in vals[1:]:
    if val == None:
      output.append("n/a")
      fmt.append("{:>13s}  ")
      continue
    if typ == 'intrinsic':
      err = val - vals[0]
      if abs(err) < tol:
        fmt.append("{:13.7f}  ")
      else:
        fmt.append("\033[1m{:13.7f}\033[0m  ")
    elif typ == 'extrinsic':
      if (vals[0] == 0 and val == 0):
        err = 0.
      elif vals[0] == 0:
        err = 1000000000.
      else:
        err = abs((val - vals[0])/vals[0]) 
      if abs(err) < tol:
        fmt.append("{:13.6%}  ")
      else:
        fmt.append("\033[1m{:13.6%}\033[0m  ")
    elif typ == 'display':
      err = val
      fmt.append("{:13.7f}  ")
    elif typ == 'ratio':
      if (vals[0] == 0 and val == 0):
        err = 1.
      else:
        err = val / vals[0]
      fmt.append("{:13.6%}  ")
    output.append(err)

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
        if len(pattern[out_name]) == 3:
          scale = float(pattern[out_name][2])
        else:
          scale = 1.
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
                val = scale*float(tmp)
              except ValueError:
                val = 0
        break
    vals.append(val)
  compare_vals(vals, name, tol, typ)
  return vals

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


def get_forces(runs):
  all_forces = []
  for run in runs:
    with open(run + "/run0.out") as f:
      f.seek(0)
      lines = f.readlines()
      for i in range(len(lines)):
        if "Forces acting on atoms" in lines[i]:
          break
      i = i + 1
      while not "atom" in lines[i]:
        i = i + 1
      forces = []
      while "atom" in lines[i]:
        toks = lines[i].split()
        forces.append(float(toks[6]))
        forces.append(float(toks[7]))
        forces.append(float(toks[8]))
        i = i + 1
    all_forces.append(forces)
  return all_forces

def rmse(data):
  rmses = []
  sqmag = [x*x for x in data[0]]
  mean = sum(sqmag) / len(sqmag)
  rmses.append(sqrt(mean))

  for i in range(1,len(data)):
    diff = [data[i][j] - data[0][j] for j in range(len(data[0]))]
    sqdiff = [x*x for x in diff]
    mean = sum(sqdiff) / len(sqdiff)
    rmses.append(sqrt(mean))
  return rmses 

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

def run_test(inputs, exe, testdir, opts):
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
  if opts.ref in runs:
    runs.remove(opts.ref)
    copytree('./'+opts.ref, './old')
    runs.insert(0, 'old')
    runs.insert(1, opts.ref)

  # mask the runs
  if opts.runs != None:
    do_runs = opts.runs.split(":")
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
      system('vasp-run.py -n ' + str(opts.nproc) + ' -e ../bin/ > /dev/null')
    else:
      system('qe-run.py ' + inputs + '.'+run+'.in --nproc ' + str(opts.nproc)  + ' --npot ' + str(opts.npot) + ' --npool ' + str(opts.npool) + ' -e ../bin/ > /dev/null')
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
                    "nb":     opts.nb,
                    "etot":   0.005,
                    "efermi": 0.1,
                    "force":  0.01,
                    "stress": 1.,
                    "bands":  0.01,
                    "atoms": 1.
                   } 
  config = dict(list(default_config.items()) + list(loaded_config.items()))
  total_energy = {"OUTCAR": ["free  energy   TOTEN", 4, 1./config["atoms"]], "run0.out": ["!    total energy", 4, 13.6/config["atoms"]]}
  fermi_energy = {"OUTCAR": ["E-fermi", 2], "run0.out": ["the Fermi energy", 4]}
  total_force = {"run0.out":  ["Total force =", 3]}
  pressure =     {"OUTCAR": ["external pressure", 3], "run0.out": ["total   stress", 5]} 
  stressXX = {"OUTCAR": ["  Total      ", 1]}
  stressYY = {"OUTCAR": ["  Total      ", 2]}
  stressZZ = {"OUTCAR": ["  Total      ", 3]}
  stressXY = {"OUTCAR": ["  Total      ", 4]}
  stressXZ = {"OUTCAR": ["  Total      ", 6]}
  stressYZ = {"OUTCAR": ["  Total      ", 5]}

  # Compare some fields
  disp_output(runs, total_energy, "Total Energy", 'intrinsic', config['etot'])
  ef = disp_output(runs, fermi_energy, "Fermi Energy", 'intrinsic', config['efermi'])
  #disp_output(runs, total_force, "Total Force", 'extrinsic', config['force'])
  disp_output(runs, pressure, "Pressure", 'intrinsic', config['stress'])
  disp_output(runs, stressXX, "Stress (XX)", 'intrinsic', config['stress'])
  disp_output(runs, stressYY, "Stress (YY)", 'intrinsic', config['stress'])
  disp_output(runs, stressZZ, "Stress (ZZ)", 'intrinsic', config['stress'])
  disp_output(runs, stressXY, "Stress (XY)", 'intrinsic', config['stress'])
  disp_output(runs, stressXZ, "Stress (XZ)", 'intrinsic', config['stress'])
  disp_output(runs, stressYZ, "Stress (YZ)", 'intrinsic', config['stress'])
  test_eigvals(runs, ef, config["bands"], config["nb"])
#  compare_vals(rmse(get_forces(runs)), "RMSE Force")

  print("------------------------------------------------------")    

  gvecs = {"OUTCAR": ["total plane-waves  NPLWV", 4], "run0.out": ["Kohn-Sham Wavefunctions", 5]}
  ecut = {"run0.out": ["kinetic-energy cutoff", 3], "OUTCAR": ["ENCUT  =  ", 4]}
  electroncs = {"run0.out": ["number of electrons       =", 4], "OUTCAR": ["NELECT =  ", 2]}
  iters = {"run0.out": ["convergence has been achieved", 5]}
  nks = {"run0.out": ["number of k points=", 4]}
  conv_thr = {"run0.out": ["convergence threshold     =", 3]}
  disp_output(runs, ecut,     "Energy Cutoff")
  disp_output(runs, gvecs,    "# Gvectors")
  disp_output(runs, iters,    "# Iters")
  disp_output(runs, nks,      "# Kpoints")
  disp_output(runs, conv_thr, "Conv Thr")


  print("------------------------------------------------------")    
  scf_time = {"OUTCAR": ["LOOP+:", 3], "run0.out": ["electrons    :", 4]}
  nscf_time = {"run1.out": ["electrons    :", 4]}
  forces_time = {"run0.out": ["forces       :", 4]}
  stresses_time = {"run0.out": ["stress       :", 4]}
  fft_time   = {"run.err": ["  FFT", 1]}
  proj_time  = {"run.err": ["  PROJ", 1]}
  vnl_time   = {"run.err": ["  VNL", 1]}
  kv_time    = {"run.err": ["  KV", 1]}
  diag_time  = {"run.err": ["  DIAG", 1]}
  FandS_time = {"run.err": ["  FandS", 1]}
  rho_time   = {"run.err": ["  RHO", 1]}

  frho_time   = {"run.err": ["  fRHO", 1]}
  fvloc_time   = {"run.err": ["  fVloc", 1]}
  fact_time   = {"run.err": ["  fAct", 1]}
  fvnl_time   = {"run.err": ["  fVnl", 1]}
  else_time  = {"run.err": ["  Else", 1]}
  disp_output(runs, scf_time, "SCF Time", 'ratio')
  disp_output(runs, nscf_time, "NSCF Time", 'ratio')
  disp_output(runs, forces_time, "Force Time", 'ratio')
  disp_output(runs, stresses_time, "Stress Time", 'ratio')

  disp_output(runs, frho_time, "fRho Time", 'display')
  disp_output(runs, fvloc_time, "fVloc Time", 'display')
  disp_output(runs, fvnl_time, "fVnl Time", 'display')
  disp_output(runs, fact_time, "fAct Time", 'display')

  disp_output(runs, fft_time, "FFT Time", 'display')
  disp_output(runs, proj_time, "PROJ Time", 'display')
  disp_output(runs, vnl_time, "VNL Time", 'display')
  disp_output(runs, kv_time, "KV Time", 'display')
  disp_output(runs, diag_time, "DIAG Time", 'display')
  disp_output(runs, FandS_time, "FandS Time", 'display')
  disp_output(runs, rho_time, "RHO Time", 'display')
  disp_output(runs, else_time, "Else Time", 'display')

  chdir(cwd)
  return 1
'''
  ef = disp_output(runs, "the Fermi energy", 4, "Fermi Energy", config['efermi'])
  disp_output(runs, "convergence has been achieved in", 5, "N Iter", 1)
  disp_output(runs, "Total force", 3, "Total Force", config['force'], intrinsic=False)
  disp_output(runs, "temperature", 2, "Temperature", 0.01)

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
parser.add_option("--ref", dest="ref", default="ref",
                  help="Use REF as the reference for comparison")
parser.add_option("-t", "--test", dest="test", default=None,
                  help="Directory containing the test")
parser.add_option("-s", "--set", dest="test_set", default=None,
                  help="File specifying set of tests to run") 
parser.add_option("--nproc", type=int, dest="nproc", default=1,
                  help="Number of MPI PEs")
parser.add_option("--npot", type=int, dest="npot", default=-1,
                  help="Number of pots (srb pools)")
parser.add_option("--npool", type=int, dest="npool", default=1,
                  help="Number of pools")
#parser.add_option("-f", "--force", action="store_true", dest="force", default=False,
#                  help="Force the reference run to occur")
#parser.add_option("-i", "--ipm", action="store_true", dest="ipm", default=False,
#                  help="Use IPM to analyze runs")
parser.add_option("-b", "--nband", type=int, dest="nb", default=-1,
                  help="How many bands to compare")


(opts, args) = parser.parse_args()

if (opts.test == None and opts.test_set == None) or opts.exe == None:
  parser.error("Must specify at least an executable and test system")

if opts.npot == -1:
  opts.npot = opts.nproc

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

# copy testcase into new directory
print("======================================================")    
for test in test_set:

  # Generate run directories 
  testdir = './internal/' + ''.join(random.choice(string.ascii_uppercase + string.digits) for x in range(4))

  print("Test Name: " + test)
  print("Run on: " + time.strftime("%Y.%m.%d") + "  In: " + testdir)


  # run first test
  time1 = run_test(test, opts.exe, testdir, opts)

  print("======================================================")    

