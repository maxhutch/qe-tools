
from files.pwscf import PwOut

def test_load_forces():
    foo = PwOut()
    foo.load("run0.out")
    assert len(foo.forces) != 0

