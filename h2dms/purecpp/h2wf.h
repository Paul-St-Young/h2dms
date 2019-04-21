#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H
#include "common.h"
class WaveFunction
{ // use constant wavefunction as base class
public:
  WaveFunction(){}
  virtual double lnwf(const Matrix &pos){return 0.0;};
  virtual Vector grad_lnwf(const Matrix &pos, const int i)
  {
    Vector vec(ndim);
    vec = Vector::Zero(ndim);
    return vec;
  };
  virtual double lap_lnwf(const Matrix &pos, const int i){return 0.0;};
};
#endif
