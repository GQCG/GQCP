#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "Atom.hpp"

namespace py = pybind11;

PYBIND11_MODULE(gqcpy, m) {
    py::class_<GQCP::Atom>(m, "Atom")
            .def(py::init<int, double, double, double>());

}
