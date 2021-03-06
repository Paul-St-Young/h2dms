#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
namespace py=pybind11;

#include "h2hp.h"
#include "h2ham.h"
#include "j2pade.h"
#include "prodwf.h"

PYBIND11_MODULE(h2hp, m)
{
  py::class_<H2HartreeProduct>(m, "H2HartreeProduct")
    .def(py::init<double, double>())
    .def_readonly("ions", &H2HartreeProduct::ions)
    .def("lnwf", &H2HartreeProduct::lnwf)
    .def("grad_lnwf", &H2HartreeProduct::grad_lnwf)
    .def("lap_lnwf", &H2HartreeProduct::lap_lnwf)
    ;
  py::class_<PadePairJastrow>(m, "PadePairJastrow")
    .def(py::init<double, double, double>())
    .def_readwrite("a", &PadePairJastrow::a)
    .def_readwrite("b", &PadePairJastrow::b)
    .def_readwrite("c", &PadePairJastrow::c)
    .def("lnwf", &PadePairJastrow::lnwf)
    .def("grad_lnwf", &PadePairJastrow::grad_lnwf)
    .def("lap_lnwf", &PadePairJastrow::lap_lnwf)
    ;
  py::class_<ProductWaveFunction>(m, "ProductWaveFunction")
    .def(py::init(&ProductWaveFunction::create_sj))
    .def("lnwf", &ProductWaveFunction::lnwf)
    .def("grad_lnwf", &ProductWaveFunction::grad_lnwf)
    .def("lap_lnwf", &ProductWaveFunction::lap_lnwf)
    ;
  py::class_<H2Hamiltonian>(m, "H2Hamiltonian")
    .def(py::init<const Matrix&, const ProductWaveFunction&>())
    .def("kinetic", &H2Hamiltonian::kinetic)
    .def("ii", &H2Hamiltonian::ii)
    .def("ei", &H2Hamiltonian::ei)
    .def("ee", &H2Hamiltonian::ee)
    .def("potential", &H2Hamiltonian::potential)
    .def("local", &H2Hamiltonian::local)
    ;
}
