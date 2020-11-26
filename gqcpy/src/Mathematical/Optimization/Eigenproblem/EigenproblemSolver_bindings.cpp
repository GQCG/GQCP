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

#include "Mathematical/Optimization/Eigenproblem/Davidson/DavidsonSolver.hpp"
#include "Mathematical/Optimization/Eigenproblem/EigenproblemSolver.hpp"

#include <pybind11/pybind11.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


void bindEigenproblemSolver(py::module& module) {

    auto module_eigenproblem_solver = module.def_submodule("EigenproblemSolver");

    module_eigenproblem_solver.def("Dense",
                                   &EigenproblemSolver::Dense,
                                   "Return an algorithm that can diagonalize a dense matrix.");

    module_eigenproblem_solver.def("Davidson",
                                   &EigenproblemSolver::Davidson,
                                   py::arg("number_of_requested_eigenpairs") = 1,
                                   py::arg("maximum_subspace_dimension") = 15,
                                   py::arg("convergence_threshold") = 1.0e-08,
                                   py::arg("correction_threshold") = 1.0e-12,
                                   py::arg("maximum_number_of_iterations") = 128,
                                   py::arg("inclusion_threshold") = 1.0e-03);
}


}  // namespace gqcpy
