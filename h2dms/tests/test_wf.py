#!/usr/bin/env python
import numpy as np

import sys
sys.path.insert(0, '..')
from pybind.h2hp import H2HartreeProduct

# test system
hp = H2HartreeProduct(1.4, 1.0)
pos = np.array([
  [-0.5, 0, 0],
  [ 0.5, 0, 0]
])
dx = 1e-3

def test_lnwf():
  assert np.isclose(-.4, hp.lnwf(pos))

def test_grad(atol=1e-6):
  for i in range(2):
    dlnwf0 = hp.grad_lnwf(pos, i)
    dlnwf = np.zeros(3)
    for idim in range(3):
      pos[i, idim] += dx
      lnwfp = hp.lnwf(pos)
      pos[i, idim] -= 2*dx
      lnwfm = hp.lnwf(pos)
      pos[i, idim] += dx
      dlnwf[idim] = (lnwfp-lnwfm)/(2*dx)
    assert np.allclose(dlnwf, dlnwf0, atol=atol)

def test_lap(atol=1e-4):
  lnwf0 = hp.lnwf(pos)
  for i in range(2):
    d2lnwf0 = hp.lap_lnwf(pos, i)
    d2lnwf = 0.0
    for idim in range(3):
      pos[i, idim] += dx
      lnwfp = hp.lnwf(pos)
      pos[i, idim] -= 2*dx
      lnwfm = hp.lnwf(pos)
      pos[i, idim] += dx
      d2lnwf += (lnwfp+lnwfm-2*lnwf0)/(dx**2)
    assert np.allclose(d2lnwf, d2lnwf0, atol=atol)

if __name__ == '__main__':
  test_lap()
# end __main__
