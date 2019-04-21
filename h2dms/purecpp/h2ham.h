#ifndef H2_HAMILTONIAN_H
#define H2_HAMILTONIAN_H
#include "common.h"
#include "h2wf.h"
class H2Hamiltonian
{
public:
  static constexpr double hbs2m=0.5, e2=1.0;  // Hartree atomic units
  H2Hamiltonian(const Matrix& ions, const WaveFunction& wf);
  double kinetic(const Matrix& pos)const;
  double ei(const Matrix& pos)const;
  double ee(const Matrix& pos)const;
  double potential(const Matrix& pos)const;
  double local(const Matrix& pos)const;
private:
  const Matrix& _ions;
  const WaveFunction& _wf;
  double _ii_pot;
};
#endif
