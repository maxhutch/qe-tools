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

# define test output
def compare_vals(vals, name, tol = 0.001, intrinsic=True):
  output = []
  if vals[0] == None:
    return
  output.append(vals[0])
  for val in vals[1:]:
    if not intrinsic:
      if (vals[0] == 0 and val == 0):
        err = 0.
      elif vals[0] == 0:
        err = 1000000000.
      else:
        err = abs((val - vals[0])/vals[0])
    else:
      err = val - vals[0]
    output.append(err)
  fmt = "{:13s}  "+"{:6.2f}  "*len(output)
  print(fmt.format(name,*output))

# define a value tester
def test_output(runs, pattern, position, name, tol, intrinsic=True):
  vals = []
  for run in runs:
    with open(run + "/run0.out", 'r') as f:
      f.seek(0)
      val = None
      for line in f:
        if pattern in line:
          toks = line.split()
          if len(toks) <= position:
            continue
          tmp = toks[position]
          if tmp[-1] == ':':
            tmp = tmp[0:-1]
          try: 
            val = float(tmp)
          except ValueError:
            val = 0
    vals.append(val)
  compare_vals(vals, name, tol, intrinsic)
  return vals

def diff_eigenvalues(out_name, ref_name, efo = 0., efr = 0., nb = -1):
  ovals = numpy.loadtxt(out_name)
  rvals = numpy.loadtxt(ref_name)

  nk = min(ovals.shape[0], rvals.shape[0]) 
  if nb < 0:
    nb = min(ovals.shape[1], rvals.shape[1]) - 4

  diff = ovals[0:nk, 4:nb+4] - rvals[0:nk, 4:nb+4] - efo + efr
  chi2 = (diff * diff).sum(axis=1) / nb
  chit = numpy.sqrt(chi2.sum() / nk)

  return chit

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

'''
# define a comparator
def compare_output(pattern, position, name, out, ref):
  out.seek(0)
  val_o = 0.
  for line in out:
    if pattern in line:
      toks = line.split()
      tmp = toks[position]
      if tmp[-1] == ':':
        tmp = tmp[0:-1]
      val_o = val_o + float(tmp)
  ref.seek(0)
  val_r = 0.
  for line in ref:
    if pattern in line:
      toks = line.split()
      tmp = toks[position]
      if tmp[-1] == ':':
        tmp = tmp[0:-1]
      val_r = val_r + float(tmp)
  print("%5.2fx  %-15s % 13.4f vs % 13.4f" % ((val_r+.0000001)/(val_o+.0000001), name, val_o, val_r))

def get_tag(f, pattern, pos):
  f.seek(0)
  for line in f:
    if pattern in line:
      toks = line.split()
      ans = toks[pos]
      return ans
'''


def run_test(inputs, exe, testdir, np = 1, ipm = False, force = False, nb = -1):
  # setup step
  print("------------------------------------------------------")    

  # Copy test into internal directory
  rmtree(testdir, ignore_errors = True)
  copytree(inputs, testdir)
  mkdir(testdir + '/bin')
  copy2(exe + '/pw.x', testdir+'/bin')
  cwd = getcwd()
  chdir(testdir)

#  print(' NP = ' + str(np))

  # Setup list of runs
  runs = next(walk('.'))[1];  
  runs = [x for x in runs if x != 'pseudo' and x != 'bin']
  # Make a copy of the 'ref' run for comparison
  if 'ref' in runs:
    runs.remove('ref')
    copytree('./ref', './old')
    runs.insert(0, 'old')
    runs.insert(1, 'ref')
 
  # Loop over the runs, running each and recording time 
  fmt = "{:13s}  "+"{:6s}  "*len(runs)
  print(fmt.format("runs:", *runs))
  print("{:13s}  ".format("time:"), end=""); stdout.flush()
  for run in runs:
    if run == 'old':
      with open('./old/time', 'r') as f:
        print("%5.2f\t" % float(f.readline()), end=""); stdout.flush()
      continue
    chdir(run)
    start = time.time()  
    system('qe-run.py ' + inputs + '.'+run+'.in -n ' + str(np) + ' -e ../bin/ > /dev/null')
    time_test = time.time() - start
    system('echo "'+str(time_test)+'" > time')
    print("%5.2f\t" % time_test , end=""); stdout.flush()
    chdir('..')
  print("")

  # This is the old test/ref comparison.  will be replaced
  print("------------------------------------------------------")    
#  with open('./test/run0.out') as test:
#   with open('./ref/run0.out') as ref:
  test_output(runs, "!    total energy", 4, "Total Energy", 0.01)
  test_output(runs, "the Fermi energy", 4, "Fermi Energy", 0.01)
  test_output(runs, "convergence has been achieved in", 5, "N Iter", 1)
  test_output(runs, "Total force", 3, "Total Force", 0.01)
  test_output(runs, "total   stress", 5, "Total Stress", 1.)
  test_output(runs, "temperature", 2, "Temperature", 0.01)
#  compare_bands = exists('./ref/bands.dat')

#  if compare_bands:
#    chi_o = diff_eigenvalues('./test/bands.dat', './ref/bands.dat', efo, efr, nb)
#    compare_vals('RMS Error', chi_o, 0.000, 0.005)

  chdir(cwd)
  return 1

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
parser.add_option("-r", "--reference", dest="exe2", default=None,
                  help="Reference executables for comparison testing")
parser.add_option("-t", "--test", dest="test", default=None,
                  help="Directory containing the test")
parser.add_option("-s", "--set", dest="test_set", default=None,
                  help="File specifying set of tests to run") 
parser.add_option("-n", "--num_pe", type=int, dest="np", default=1,
                  help="Number of MPI PEs")
parser.add_option("-f", "--force", action="store_true", dest="force", default=False,
                  help="Force the reference run to occur")
parser.add_option("-i", "--ipm", action="store_true", dest="ipm", default=False,
                  help="Use IPM to analyze runs")
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
  time1 = run_test(test, opts.exe, testdir, opts.np, opts.ipm, opts.force, opts.nb)

  # if on a comparison, run reference 
  if opts.exe2 != None:
    time2 = run_test(test, opts.exe2, testdir2, opts.np, opts.ipm, opts.nb) 

  print("======================================================")    

