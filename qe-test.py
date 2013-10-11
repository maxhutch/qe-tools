#!/usr/bin/python3

from sys import argv, exit
from os import system, chdir, mkdir, getcwd
from os.path import exists
from shutil import move, copytree, copy2, rmtree
import time
from optparse import OptionParser
import random
import string 
import numpy

# define test output
def compare_vals(name, oval, rval, tol = 0.001, intrinsic=True):
  if not intrinsic:
    if (rval == 0 and oval == 0):
      err = 0.
    elif rval == 0:
      err = 1000000000.
    else:
      err = abs((oval - rval)/rval)
  else:
    err = oval - rval
  if abs(err) < tol:
    print(" passed  %-14s % 10.6e vs % 10.6e  (% 10.8f) " % (name, oval, rval, err))
  else:
    print(" FAILED  %-14s % 10.6e vs % 10.6e  (% 10.8f) " % (name, oval, rval, err))

# define a value tester
def test_output(pattern, position, name, out, ref, tol, intrinsic=True):
  out.seek(0)
  val_r = None
  for line in out:
    if pattern in line:
      toks = line.split()
      if len(toks) <= position:
        continue
      tmp = toks[position]
      if tmp[-1] == ':':
        tmp = tmp[0:-1]
      try: 
        val_o = float(tmp)
      except ValueError:
        val_o = 0
  ref.seek(0)
  for line in ref:
    if pattern in line:
      toks = line.split()
      if len(toks) <= position:
        continue
      tmp = toks[position]
      if tmp[-1] == ':':
        tmp = tmp[0:-1]
      try:
        val_r = float(tmp)
      except ValueError:
        val_r = 0
  if val_r == None:
    return (0.,0.)
  compare_vals(name, val_o, val_r, tol, intrinsic)
  return (val_o, val_r)

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
  # Do we need to do a reference run?
  rerun = not exists(inputs + '/ref') or force
  if not rerun:
    system('diff -w ' + inputs + '/ref/*.ref.in ' + inputs + '/*.ref.in > tmp')
    with open('./tmp', 'r') as f:
      if len(f.readlines()) != 0:
        rerun = True
      else: 
        print('Using saved result')
    system('rm tmp')

  rmtree(testdir, ignore_errors = True)
  copytree(inputs, testdir)
  mkdir(testdir + '/bin')
  copy2(exe + '/pw.x', testdir)
  copy2(exe + '/pp.x', testdir)
  copy2(exe + '/plotrho.x', testdir)
#  copy2(exe + '../scripts/run.py', testdir)
  cwd = getcwd()
  chdir(testdir)

  print(' NP = ' + str(np))

  start = time.time()
  system('qe-run.py ' + inputs + '.test.in -n ' + str(np) + ' > /dev/null')
  time_test = time.time() - start

  compare_bands = exists('./bands.dat')

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

  print("------------------------------------------------------")    
  with open('./run0.test.out') as test:
   with open('./run0.ref.out') as ref:
     test_output("!    total energy", 4, "Total Energy", test, ref, 0.01)
     (efo, efr) = test_output("the Fermi energy", 4, "Fermi Energy", test, ref, 0.01)
     test_output("convergence has been achieved in", 5, "N Iter", test, ref, 1)
     test_output("Total force", 3, "Total Force", test, ref, 0.01)
     test_output("total   stress", 5, "Total Stress", test, ref, 1.)
     test_output("temperature", 2, "Temperature", test, ref, 0.01)

  if compare_bands:
    chi_o = diff_eigenvalues('./bands.test.dat', './bands.ref.dat', efo, efr, nb)
    compare_vals('RMS Error', chi_o, 0.000, 0.005)

  chdir(cwd)

  return 1

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

