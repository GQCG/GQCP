#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "Atom.hpp"
#include "HamiltonianParameters/HamiltonianParameters.hpp"
#include "Molecule.hpp"
#include "typedefs.hpp"

namespace py = pybind11;

PYBIND11_MODULE(gqcpy, m) {
    py::class_<GQCP::Atom>(m, "Atom")
            .def(py::init<int, double, double, double>())
            ;
    py::class_<GQCP::Molecule>(m, "Molecule")
            .def(py::init<std::vector<GQCP::Atom>, int>())
            .def("read_xyz", &GQCP::Molecule::Readxyz)
            .def("get_N", &GQCP::Molecule::get_N)
            .def("get_atoms", &GQCP::Molecule::get_atoms)
            .def("number_of_atoms", &GQCP::Molecule::numberOfAtoms)
            ;
}