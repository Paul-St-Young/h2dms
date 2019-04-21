import numpy as np

import sys
sys.path.insert(0, '..')
from pybind.h2hp import H2HartreeProduct

def test_grad():
  hp = H2HartreeProduct(1.4, 1.0)
  pos = np.array([
    [-0.5, 0, 0],
    [ 0.5, 0, 0]
  ])
  dx = 1e-3
  for i in range(2):
    dlnwf0 = hp.diff_lnwf(pos, i)
    dlnwf = np.zeros(3)
    for idim in range(3):
      pos[i, idim] += dx
      lnwfp = hp.lnwf(pos)
      pos[i, idim] -= 2*dx
      lnwfm = hp.lnwf(pos)
      dlnwf[idim] = (lnwfp-lnwfm)/(2*dx)
    assert np.allclose(dlnwf, dlnwf0)
