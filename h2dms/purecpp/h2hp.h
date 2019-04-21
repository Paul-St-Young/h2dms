#ifndef H2_HARTREE_PRODUCT_H
#define H2_HARTREE_PRODUCT_H
#include <Eigen/Dense>


class H2HartreeProduct
{
public:
  typedef Eigen::VectorXd Vector;
  typedef Eigen::MatrixXd Matrix;
  Matrix ions;
  static const int natom=2, ndim=3;
  H2HartreeProduct(double rbond, double alpha);
  // wf value
  double lnwf(const Matrix& pos);
  // wf ratio
  Vector grad_lnwf(const Matrix& pos, const int i);
  double lap_lnwf(const Matrix& pos, const int i);
  // getters
  double get_rbond(){return _rbond;}
  // setters
  void set_rbond(double rbond){_rbond=rbond;}
private:
  double _rbond;
  double _alpha;
};
#endif
