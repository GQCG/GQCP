// This file is part of GQCG-GQCP.
//
// Copyright (C) 2017-2020  the GQCG developers
//
// GQCG-GQCP is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// GQCG-GQCP is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-GQCP.  If not, see <http://www.gnu.org/licenses/>.

#include "Mathematical/Algorithm/IterativeAlgorithm.hpp"
#include "QCMethod/Geminals/AP1roG.hpp"

#include <pybind11/pybind11.h>


namespace py = pybind11;


namespace gqcpy {


void bindQCMethodAP1roG(py::module& module) {
    py::class_<GQCP::QCMethod::AP1roG>(module, "AP1roG", "The AP1roG quantum chemical method.")

        // CONSTRUCTORS

        .def(py::init<GQCP::RSQHamiltonian<double>, size_t>(),
             py::arg("sq_hamiltonian"),
             py::arg("N_P"))


        // PUBLIC METHODS

        .def(
            "optimize",
            [](const GQCP::QCMethod::AP1roG& qc_method, GQCP::IterativeAlgorithm<GQCP::NonLinearEquationEnvironment<double>>& solver, GQCP::NonLinearEquationEnvironment<double>& environment) {
                return qc_method.optimize(solver, environment);
            },
            "Optimize the AP1roG wave function model.");
}


}  // namespace gqcpy
