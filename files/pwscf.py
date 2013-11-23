import numpy as np

class PwOut(object):
  def __init__(self, fname = None):
    self.forces = np.zeros((0,3))
    if fname != None:
      self.load(fname)

  def load(self,fname):
    with open(fname, 'r') as f:
      lines = f.readlines()

    force_line = 0
    for i in range(len(lines)):
      if "Forces acting on atoms (Ry/au):" in lines[i]:
        force_line = i + 4
        break
    while True: 
      toks = lines[force_line].split()
      if len(toks) != 9:
        break
      self.forces = np.append(self.forces, [[float(toks[6]), float(toks[7]), float(toks[8])]], axis=0)
      force_line = force_line + 1
    return 0

