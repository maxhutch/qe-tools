qe-tools
========

qe-tools is a collection of tools for interacting with Quantum Espresso.
There are two types of tools so far:

 * Tools which will run QE for you
 * Tools to process QE outputs

List of tools
-------------

### [run.py](run.py)

This script parses and runs a special input file.  It reads the file sequentially,
adding or overwriting namelist options and keeping only the most recently read
card of each type.  When it hits a "@run" directive, it dumps state as a valid
pwscf input file and runs pwscf.  If there are multiple runs in a single file, 
their stdout/stderr are directed to sequentially numbered run.out files.

### [compare_bands.py](compare_bands.py)

This script takes a reference and test band-structure file, and optionally band
ranges and Fermi energies, computes the RMS error between the band-structures
and plots with gnuplot.

### [test.py](test.py)

This script runs validation tests.  If the test is named NAME, it expects
a directory NAME containing NAME.ref.in and NAME.test.in files, which it 
will run and compare on the basis of energy, Fermi energy, force, stress,
and RMS band-stucture error (all if applicable).  The input files are of the
same format as compare\_bands.py.  There can optionally be a directory ref
in the same directory as the inputs which contains a previous reference output.
If the ref directory is present, its contents are reused instead of re-running
the reference.

