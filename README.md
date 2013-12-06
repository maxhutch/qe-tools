qe-tools
========

qe-tools is a collection of tools for interacting with Quantum Espresso.
There are two types of tools so far:

 * Tools which will run QE for you
 * Tools to process QE outputs

List of tools
-------------

### [qe-test.py](qe-test.py)

This script runs validation tests. If the test is named TEST, it expects
runs TEST/RUN with a qe-run input file TEST.RUN.in.  A special run is called
'ref' and should also contain the outputs of qe-run from a trusted build.
The runs are executed and outputs compared to the 'ref' run.  The 'ref' run
is also re-executed to validate the build.  Here is [an example](https://github.com/maxhutch/qe-tests/tree/master/C).

This script can also handle VASP.  The structure is the same as for QE, but
with the contents of TEST/RUN being the VASP input files, usually just 
INCAR, POSCAR, POTCAR, KPOINTS.

### [qe-run.py](qe-run.py)

This script parses and runs a special input file.  It reads the file sequentially,
adding or overwriting namelist options and keeping only the most recently read
card of each type.  When it hits a "@run" directive, it dumps state as a valid
pwscf input file and runs pwscf.  If there are multiple runs in a single file, 
their stdout/stderr are directed to sequentially numbered run.out files.

### [vasp-run.py](vasp-run.py)

This script is a thin wrapper around the command:
  
    mpirun -np N BINDIR/vasp 2> PREFIX.err | tee PREFIX.out

which pulls N from the option -n, the BINDIR from the option -e, and the PREFIX
from the -p option, as in [qe-run.py](qe-run.py).

### [compare_bands.py](compare_bands.py)

This script takes a reference and test band-structure file, and optionally band
ranges and Fermi energies, computes the RMS error between the band-structures
and plots with gnuplot.



