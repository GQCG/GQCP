
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "Atom.hpp"
#include "HamiltonianParameters/HamiltonianParameters.hpp"
#include "CISolver/CISolver.hpp"
#include "HamiltonianBuilder/Hubbard.hpp"
#include "Drivers/HubbardDriver.hpp"

namespace py = pybind11;

PYBIND11_MODULE(gqcpy, m) {
    py::class_<GQCP::Atom>(m, "Atom")
            .def(py::init<int, double, double, double>());

    py::class_<GQCP::HubbardDriver>(m, "HubbardDriver")
            .def(py::init<std::string, size_t, size_t, size_t, size_t>())
            .def("get_energies", &GQCP::HubbardDriver::get_energies);

}
