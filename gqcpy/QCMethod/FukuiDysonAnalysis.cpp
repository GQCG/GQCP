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
#include "QCMethod/FukuiDysonAnalysis.hpp"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>


namespace py = pybind11;


namespace gqcpy {


void bindFukuiDysonAnalysis(py::module& module) {
    py::class_<GQCP::QCMethod::FukuiDysonAnalysis>(module, "FukuiDysonAnalysis", "A class that solves the FCI Hamiltonian for a given molecule and performs Fukui and Dyson analysis")
        .def(py::init<const GQCP::Molecule& , const std::string&>(), py::arg("molecule"), py::arg("basis_set"))
        .def("get_dyson_coefficients", &GQCP::QCMethod::FukuiDysonAnalysis::get_dyson_coefficients)
        .def("get_fukui_matrix", &GQCP::QCMethod::FukuiDysonAnalysis::get_fukui_matrix)
        .def("get_fukui_naturals", &GQCP::QCMethod::FukuiDysonAnalysis::get_fukui_naturals)
        .def("get_fukui_vectors", &GQCP::QCMethod::FukuiDysonAnalysis::get_fukui_vectors)
        .def("get_canonical_matrix", &GQCP::QCMethod::FukuiDysonAnalysis::get_canonical_matrix);
}



}  // namespace gqcpy
