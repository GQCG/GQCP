#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <string>

#include "Atom.hpp"
#include "CISolver/CISolver.hpp"
#include "FockSpace/ProductFockSpace.hpp"
#include "HamiltonianBuilder/Hubbard.hpp"
#include "HamiltonianParameters/HamiltonianParameters.hpp"
#include "math/SquareMatrix.hpp"
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

/*    py::class_<GQCP::SquareMatrix<double>>(m, "SquareMatrix")
            .def(py::init())
            ;

    py::class_<GQCP::HoppingMatrix>(m, "HoppingMatrix")
            .def(py::init<const GQCP::SquareMatrix<double> &>())
            .def_static("from_upper_triangle", &GQCP::HoppingMatrix::FromUpperTriangle)
            ;*/

/*    py::class_<GQCP::HamiltonianParameters<double>>(m, "HamiltonianParameters")
            .def_static("hubbard", &GQCP::HamiltonianParameters<double>::Hubbard)
            ;*/

    py::class_<GQCP::ProductFockSpace>(m, "ProductFockSpace")
            .def(py::init<size_t, size_t, size_t>())
            ;
    py::class_<GQCP::Hubbard>(m, "Hubbard")
            .def(py::init<const GQCP::ProductFockSpace &>())
            ;
    py::class_<GQCP::CISolver>(m, "CISolver")
            .def(py::init<const GQCP::HamiltonianBuilder &, const GQCP::HamiltonianParameters<double>>())
            .def("solve", &GQCP::CISolver::solve)
            ;

    py::class_<GQCP::DenseSolverOptions>(m, "DenseSolverOptions");

}