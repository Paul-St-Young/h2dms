#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H
#include "common.h"
class WaveFunction
{ // use constant wavefunction as base class
  //  throw exception if base class functions are called
  //   i.e. block confusing polymorphism
public:
  WaveFunction(){}
  virtual double lnwf(const Matrix &pos)
  const{
    throw 0;
    return 0.0;
  };
  virtual Vector grad_lnwf(const Matrix &pos, const int i)
  const {
    throw 0;
    Vector vec(ndim);
    vec = Vector::Zero(ndim);
    return vec;
  };
  virtual double lap_lnwf(const Matrix &pos, const int i)
  const {
    throw 0;
    return 0.0;
  };
};
#endif
