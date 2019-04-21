#include "h2hp.h"

H2HartreeProduct::H2HartreeProduct(
    double rbond,  double alpha)
  : _rbond(rbond), _alpha(alpha)
{
  ions.resize(natom, ndim);
  ions << -rbond/2., 0, 0,
           rbond/2., 0, 0;
}

double H2HartreeProduct::lnwf(const Matrix& pos)
{
  double val = 0.0;
  for (int iatom=0; iatom<natom; iatom++)
  {
    val += (pos.row(iatom) - ions.row(iatom)).norm();
  }
  return -_alpha*val;
}

H2HartreeProduct::Vector
H2HartreeProduct::grad_lnwf(const Matrix& pos, const int i)
{
  Vector rvec(ndim);
  rvec = pos.row(i)-ions.row(i);
  return -_alpha*rvec/rvec.norm();
}
