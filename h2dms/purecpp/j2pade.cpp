#include "j2pade.h"

PadePairJastrow::PadePairJastrow(double a_in, double b_in, double c_in)
  : a(a_in), b(b_in), c(c_in)
{
}

double PadePairJastrow::fpade(const double r)
const {
  double nume, deno;
  nume = a*r+b*r*r;
  deno = 1+c*r;
  return nume/deno;
}

double PadePairJastrow::dfpade(const double r)
const {
  double arg, nume, deno;
  arg = 1.+c*r;
  nume = a+b*r*(1+arg);
  deno = arg*arg;
  return nume/deno;
}

double PadePairJastrow::d2fpade(const double r)
const {
  double arg, nume, deno;
  arg = 1.+c*r;
  nume = 2*(b-a*c);
  deno = arg*arg*arg;
  return nume/deno;
}

double PadePairJastrow::lnwf(const Matrix& pos)
const {
  double r, ujas=0.0;
  int nptcl = pos.rows();
  for (int i=0; i<nptcl; i++)
  {
    for (int j=i+1; j<nptcl; j++)
    {
      r = (pos.row(i)-pos.row(j)).norm();
      ujas += fpade(r);
    }
  }
  return -ujas;
}

Vector PadePairJastrow::grad_lnwf(const Matrix& pos, const int i)
const {
  int nptcl = pos.rows();
  double r=0.0;
  Vector grad(ndim), rvec(ndim);
  grad = Vector::Zero(ndim);
  rvec = Vector::Zero(ndim);
  for (int j=0; j<nptcl; j++)
  {
    if (i==j) continue;
    rvec = pos.row(i)-pos.row(j);
    r = rvec.norm();
    grad += rvec/r*dfpade(r);
  }
  return -grad;
}
