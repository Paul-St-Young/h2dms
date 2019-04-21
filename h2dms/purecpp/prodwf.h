#ifndef PRODUCT_WAVEFUNCTION_H
#define PRODUCT_WAVEFUNCTION_H
#include <list>
#include "h2wf.h"
class ProductWaveFunction : public WaveFunction
{
public:
  typedef std::list<WaveFunction*> WFList;
  ProductWaveFunction(WFList wfs);
  double lnwf(const Matrix& pos) const;
  Vector grad_lnwf(const Matrix& pos, const int i) const;
  double lap_lnwf(const Matrix& pos, const int i) const;
  // factory method
  static ProductWaveFunction create_sj(double rbond, double alpha,
    double a, double b, double c);
private:
  WFList _wfs;
};
#endif
