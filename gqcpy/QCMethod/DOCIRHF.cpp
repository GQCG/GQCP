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
#include "QCMethod/DOCIRHF.hpp"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>


namespace py = pybind11;


namespace gqcpy {


void bindQCMethodDOCIRHF(py::module& module) {
    py::class_<GQCP::QCMethod::DOCIRHF>(module, "DOCIRHF", "Finds the RHF solution and solves the DOCI Hamiltonian in that RHF basis")
        .def(py::init<const std::string&, const std::string&, const bool>(), py::arg("xyz_filename"), py::arg("basis_set"), py::arg("use_davidson"))
        .def(py::init<const GQCP::Molecule&, const std::string&, const bool>(), py::arg("molecule"), py::arg("basis_set"), py::arg("use_davidson"))
        .def("solve", &GQCP::QCMethod::DOCIRHF::solve, "Solve the eigenvalue equations such that the lowest energy and corresponding eigenvector becomes available. ")
        .def("get_energy", &GQCP::QCMethod::DOCIRHF::energy, "Get the lowest energy.")
        .def("get_rhf_energy", &GQCP::QCMethod::DOCIRHF::energy_rhf, "Get the RHF energy.")
        .def("get_transformation_matrix", &GQCP::QCMethod::DOCIRHF::transformationMatrix, "Get the total transformation matrix.");
}



}  // namespace gqcpy
