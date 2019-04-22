#!/usr/bin/env python
import numpy as np

import sys
sys.path.insert(0, '..')
from pybind.h2hp import H2HartreeProduct, H2Hamiltonian
from pybind.h2hp import ProductWaveFunction

# test system
hp = H2HartreeProduct(1.4, 1.0)
pos = np.array([
  [-0.5, 0, 0],
  [ 0.5, 0, 0]
])
dx = 1e-3

def test_kinetic():
  wf = ProductWaveFunction(1.4, 1.0, 0, 0, 0)
  ham = H2Hamiltonian(hp.ions, wf)
  assert np.isclose(9, ham.kinetic(pos))

def test_potential():
  wf = ProductWaveFunction(1.4, 1.0, 0, 0, 0)
  ham = H2Hamiltonian(hp.ions, wf)
  ii = 1./1.4
  ei = -2*(1./0.2+1./1.2)
  ee = 1.0
  assert np.isclose(ii, ham.ii(pos))
  assert np.isclose(ei, ham.ei(pos))
  assert np.isclose(ee, ham.ee(pos))
  pot = ii+ei+ei
  assert np.isclose(ii+ei+ee, ham.potential(pos))

if __name__ == '__main__':
  test_kinetic()
  test_potential()
# end __main__
