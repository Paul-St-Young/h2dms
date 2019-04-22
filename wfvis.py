#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

from h2dms.pybind.h2hp import ProductWaveFunction
from h2dms.pybind.h2hp import H2Hamiltonian

if __name__ == '__main__':
  # test system
  rbond = 1.4
  alpha = 1.0
  a = -0.5
  b = 0.0
  c = 0.0
  ions = np.array([
    [-rbond/2., 0, 0],
    [ rbond/2., 0, 0]
  ])
  wf = ProductWaveFunction(rbond, alpha, a, b, c)
  pos = np.array([
    [-0.5, 0, 0],
    [ 0.5, 0, 0]
  ])
  ham = H2Hamiltonian(ions, wf)
  # plot local energy
  myx = np.linspace(-2.0, 2.0, 64)
  yl = []
  for x in myx:
    pos[1, 0] = x
    yl.append(ham.local(pos))
  myy = np.array(yl)
  fig, ax = plt.subplots(1, 1)
  ax.plot(myx, myy)
  plt.show()
# end __main__
