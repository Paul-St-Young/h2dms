#include "h2ham.h"

H2Hamiltonian::H2Hamiltonian(const Matrix& ions, const WaveFunction& wf)
  : _ions(ions), _wf(wf), _ii_pot(0)
{
  double r = (ions.row(0)-ions.row(1)).norm();
  _ii_pot = e2/r;
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

double H2Hamiltonian::ei(const Matrix& pos)
const {
  double pot=0;
  for (int i=0; i<nelec; i++)
  {
    for (int iatom=0; iatom<natom; iatom++)
    {
      double r = (pos.row(i) - _ions.row(iatom)).norm();
      pot -= e2/r;
    }
  }
  return pot;
}

double H2Hamiltonian::ee(const Matrix& pos)
const {
  double pot=0;
  for (int i=0; i<nelec; i++)
  {
    for (int j=i+1; j<nelec; j++)
    {
      double r = (pos.row(i) - pos.row(j)).norm();
      pot += e2/r;
    }
  }
  return pot;
}

double H2Hamiltonian::potential(const Matrix& pos)
const {
  return ei(pos) + ee(pos) + _ii_pot;
}

double H2Hamiltonian::local(const Matrix& pos)
const {
  return kinetic(pos) + potential(pos);
}
