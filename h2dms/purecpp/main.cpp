#include <iostream>
#include "common.h"
#include "h2hp.h"

using namespace std;

int main(int argc, char* argv[])
{
  double rbond = 1.4; // bohr
  double alpha = 1.0; // Slater orbital exponent
  H2HartreeProduct hp(rbond, alpha);
  Matrix pos(natom, ndim);
  pos << -0.5, 0, 0,
          0.5, 0, 0;
  cout << hp.grad_lnwf(pos, 0) << endl;
  cout << hp.grad_lnwf(pos, 1) << endl;
  cout << hp.lap_lnwf(pos, 0) << endl;
  cout << hp.lap_lnwf(pos, 1) << endl;
  return 0;
}
