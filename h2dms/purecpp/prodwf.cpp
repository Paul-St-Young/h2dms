#include "prodwf.h"
// --- should move include specific wf and factory to a factory class
#include "h2hp.h"
#include "j2pade.h"

ProductWaveFunction ProductWaveFunction::create_sj(
  double rbond, double alpha,
  double a, double b, double c)
{
  H2HartreeProduct *hp = new H2HartreeProduct(rbond, alpha);
  PadePairJastrow *j2 = new PadePairJastrow(a, b, c);
  ProductWaveFunction::WFList wfs = {hp, j2};
  return ProductWaveFunction(wfs);
}
// should move include specific wf and factory to a factory class ---

ProductWaveFunction::ProductWaveFunction(
  ProductWaveFunction::WFList wfs
)
  : _wfs(wfs)
{
}

double ProductWaveFunction::lnwf(const Matrix& pos)
const {
  double val=0.0;
  for (auto it=_wfs.begin(); it!=_wfs.end(); it++)
  {
    val += (*it)->lnwf(pos);
  }
  return val;
}

Vector ProductWaveFunction::grad_lnwf(const Matrix& pos, const int i)
const {
  Vector grad(ndim);
  grad = Vector::Zero(ndim);
  for (auto it=_wfs.begin(); it!=_wfs.end(); it++)
  {
    grad += (*it)->grad_lnwf(pos, i);
  }
  return grad;
}

double ProductWaveFunction::lap_lnwf(const Matrix& pos, const int i)
const {
  double lap=0.0;
  for (auto it=_wfs.begin(); it!=_wfs.end(); it++)
  {
    lap += (*it)->lap_lnwf(pos, i);
  }
  return lap;
}
