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
#include "Mathematical/Optimization/Eigenproblem/Davidson/DavidsonSolver.hpp"
#include "Mathematical/Optimization/Eigenproblem/EigenproblemSolver.hpp"

#include <pybind11/pybind11.h>


namespace py = pybind11;


namespace gqcpy {


void bindEigenproblemSolver(py::module& module) {

    auto module_eigenproblem_solver = module.def_submodule("EigenproblemSolver");

    module_eigenproblem_solver.def("Dense",
        &GQCP::EigenproblemSolver::Dense,
        "Return an algorithm that can diagonalize a dense matrix."
    );

    module_eigenproblem_solver.def("Davidson",
        &GQCP::EigenproblemSolver::Davidson,
        py::arg("number_of_requested_eigenpairs") = 1,
        py::arg("maximum_subspace_dimension") = 15,
        py::arg("convergence_threshold") = 1.0e-08,
        py::arg("correction_threshold") = 1.0e-12,
        py::arg("maximum_number_of_iterations") = 128,
        py::arg("inclusion_threshold") = 1.0e-03
    );
}


}  // namespace gqcpy
