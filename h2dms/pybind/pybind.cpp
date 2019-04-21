#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
namespace py=pybind11;

#include "h2hp.h"

PYBIND11_MODULE(h2hp, m)
{
  py::class_<H2HartreeProduct>(m, "H2HartreeProduct")
    .def(py::init<double, double>())
    .def("lnwf", &H2HartreeProduct::lnwf)
    .def("grad_lnwf", &H2HartreeProduct::grad_lnwf)
    ;
}
