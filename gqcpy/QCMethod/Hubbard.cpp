// This file is part of GQCG-gqcp.
// 
// Copyright (C) 2017-2019  the GQCG developers
// 
// GQCG-gqcp is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// GQCG-gqcp is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-gqcp.  If not, see <http://www.gnu.org/licenses/>.
// 
#include "QCMethod/Hubbard.hpp"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>


namespace py = pybind11;


namespace gqcpy {


void bindQCMethodHubbard(py::module& module) {
    py::class_<GQCP::QCMethod::Hubbard>(module, "Hubbard", "Construct and solve the Hubbard Hamiltonian for a given upper triangular part of the (weighted) connectivity matrix.")
            .def(py::init<const std::string&, const size_t, const size_t, const size_t>(), py::arg("upper triangular part as csv"), py::arg("num_states"), py::arg("num_alpha"), py::arg("num_beta"))
            .def("solve", &GQCP::QCMethod::Hubbard::solve, "Solve the eigenvalue equations such that the energies and eigenvectors become available.")
            .def("get_energies", &GQCP::QCMethod::Hubbard::energies, "Get the set of lowest energies.")
            .def("get_one_rdms", &GQCP::QCMethod::Hubbard::oneRDMs, "Get the first order RDM corresponding to the lowest energies.");
}


}  // namespace gqcpy
