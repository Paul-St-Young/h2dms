#ifndef PADE_PAIR_JASTROW_H
#define PADE_PAIR_JASTROW_H
#include "h2wf.h"
class PadePairJastrow : public WaveFunction
{
public:
  double a, b, c;
  PadePairJastrow(double a_in, double b_in, double c_in);
  double fpade(const double r) const;
  double dfpade(const double r) const;
  double d2fpade(const double r) const;
  double lnwf(const Matrix& pos)const;
  Vector grad_lnwf(const Matrix& pos, const int i)const;
};
#endif
