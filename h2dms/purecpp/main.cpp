#include <iostream>
#include "common.h"
#include "h2hp.h"
#include "h2ham.h"
#include "j2pade.h"

using namespace std;

int main(int argc, char* argv[])
{
  double rbond = 1.4; // bohr
  double alpha = 1.0; // Slater orbital exponent
  H2HartreeProduct hp(rbond, alpha);
  Matrix pos(natom, ndim);
  pos << -0.5, 0, 0,
          0.5, 0, 0;
  //cout << hp.grad_lnwf(pos, 0).transpose() << endl;
  //cout << hp.grad_lnwf(pos, 1).transpose() << endl;
  //cout << hp.lap_lnwf(pos, 0) << endl;
  //cout << hp.lap_lnwf(pos, 1) << endl;
  Matrix ions(hp.ions);
  H2Hamiltonian ham(ions, hp);
  //cout << ham.kinetic(pos) << endl;
  //cout << ham.ii(pos) << endl;
  //cout << ham.ei(pos) << endl;
  //cout << ham.ee(pos) << endl;
  //cout << ham.potential(pos) << endl;
  //cout << ham.local(pos) << endl;
  PadePairJastrow jee(1.0, 0.1, 0.0);
  cout << jee.lnwf(pos) << endl;
  cout << jee.grad_lnwf(pos, 0).transpose() << endl;
  cout << jee.grad_lnwf(pos, 1).transpose() << endl;
  cout << jee.lap_lnwf(pos, 0) << endl;
  cout << jee.lap_lnwf(pos, 1) << endl;
  return 0;
}
