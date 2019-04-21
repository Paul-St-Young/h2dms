#include "h2ham.h"

H2Hamiltonian::H2Hamiltonian(const WaveFunction& wf)
  : _wf(wf)
{
}

double H2Hamiltonian::kinetic(const Matrix& pos)
const {
  double lap=0, kin=0;
  Vector grad(ndim);
  grad = Vector::Zero(ndim);
  for (int i=0; i<nelec; i++)
  {
    grad = _wf.grad_lnwf(pos, i);
    lap = _wf.lap_lnwf(pos, i);
    kin += hbs2m*(lap+grad.squaredNorm());
  }
  return kin;
}

double H2Hamiltonian::potential(const Matrix& pos)
const {
  return 0.0;
}

double H2Hamiltonian::local(const Matrix& pos)
const {
  return kinetic(pos) + potential(pos);
}
