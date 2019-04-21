#ifndef H2_HARTREE_PRODUCT_H
#define H2_HARTREE_PRODUCT_H
#include "h2wf.h"
class H2HartreeProduct : public WaveFunction
{
public:
  Matrix ions;
  H2HartreeProduct(double rbond, double alpha);
  // wf value
  double lnwf(const Matrix& pos) const;
  // wf ratio
  Vector grad_lnwf(const Matrix& pos, const int i) const;
  double lap_lnwf(const Matrix& pos, const int i) const;
  // getters
  double get_rbond()const{return _rbond;}
  // setters
  void set_rbond(double rbond){_rbond=rbond;}
private:
  double _rbond;
  double _alpha;
};
#endif
