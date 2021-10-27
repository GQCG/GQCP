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

#include "Mathematical/Optimization/Eigenproblem/GeneralizedEigenproblemSolver.hpp"
#include "Utilities/complex.hpp"

#include <pybind11/pybind11.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


void bindGeneralizedEigenproblemSolver(py::module& module) {

    auto module_generalized_eigenproblem_solver = module.def_submodule("GeneralizedEigenproblemSolver");

    module_generalized_eigenproblem_solver.def("Dense_d",
                                               &GeneralizedEigenproblemSolver::Dense<double>,
                                               "Return an algorithm that can diagonalize a real-valued dense matrix.");

    module_generalized_eigenproblem_solver.def("Dense_cd",
                                               &GeneralizedEigenproblemSolver::Dense<complex>,
                                               "Return an algorithm that can diagonalize a complex-valued dense matrix.");
}


}  // namespace gqcpy
