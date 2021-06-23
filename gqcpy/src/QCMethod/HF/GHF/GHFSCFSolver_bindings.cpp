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

#include "QCMethod/HF/GHF/GHFSCFSolver.hpp"
#include "gqcpy/include/interfaces.hpp"

#include <pybind11/pybind11.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


/*
 *  Add Python bindings for GHF SCF solvers.
 */
void bindGHFSCFSolvers(py::module& module) {

    // Provide bindings for real-valued GHF SCF solvers.
    py::class_<GHFSCFSolver<double>> py_GHFSCFSolver_d {module, "GHFSCFSolver_d", "A factory that can create real-valued GHF SCF solvers."};

    bindHartreeFockSCFSolverInterface(py_GHFSCFSolver_d);


    // Provide bindings for complex-valued GHF SCF solvers.
    py::class_<GHFSCFSolver<complex>> py_GHFSCFSolver_cd {module, "GHFSCFSolver_cd", "A factory that can create complex-valued GHF SCF solvers."};

    bindHartreeFockSCFSolverInterface(py_GHFSCFSolver_cd);
}


}  // namespace gqcpy
