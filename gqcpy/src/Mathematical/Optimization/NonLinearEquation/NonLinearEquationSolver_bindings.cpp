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

#include "Mathematical/Optimization/NonLinearEquation/NonLinearEquationSolver.hpp"

#include <pybind11/pybind11.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


void bindNonLinearEquationSolver(py::module& module) {

    auto module_non_linear_equation_solver = module.def_submodule("NonLinearEquationSolver");

    module_non_linear_equation_solver.def(
        "Newton",
        [](const double threshold = 1.0e-08, const size_t maximum_number_of_iterations = 128) {
            return NonLinearEquationSolver<double>::Newton(threshold, maximum_number_of_iterations);
        },
        py::arg("threshold") = 1.0e-08,
        py::arg("maximum_number_of_iterations") = 128,
        "Return an iterative algorithm that performs Newton steps to solve a system of equations.");
}


}  // namespace gqcpy
