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

#include "Mathematical/Algorithm/Algorithm.hpp"
#include "Mathematical/Algorithm/IterativeAlgorithm.hpp"
#include "QCMethod/CI/DOCINewtonOrbitalOptimizer.hpp"

#include <pybind11/pybind11.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


/**
 *  Bind a templated DOCI Newton orbital optimizer.
 * 
 *  @tparam EigenproblemSolver          the type of the eigenvalue problem solver that should be used
 * 
 *  @param module                       the pybind11 module in which these bindings should appear
 *  @param suffix                       the suffix that should appear in the Python class name in the following manner: "DOCINewtonOrbitalOptimizer_suffix"
 */
template <typename EigenproblemSolver>
void bindDOCINewtonOrbitalOptimizer(py::module& module, const std::string& suffix) {

    py::class_<DOCINewtonOrbitalOptimizer<EigenproblemSolver>>(module,
                                                               ("DOCINewtonOrbitalOptimizer_" + suffix).c_str(),
                                                               "A class that performs gradient-and-Hessian-based orbital optimization for DOCI")

        // PUBLIC METHODS

        .def(
            "optimize",
            [](DOCINewtonOrbitalOptimizer<EigenproblemSolver>& optimizer, RSpinOrbitalBasis<double, GTOShell>& spinor_basis, RSQHamiltonian<double>& sq_hamiltonian) {
                optimizer.optimize(spinor_basis, sq_hamiltonian);
            },
            py::arg("spinor_basis"),
            py::arg("sq_hamiltonian"));
}


/**
 *  Bind all kinds of DOCI Newton orbital optimizers to the given module.
 * 
 *  @param module           the module to which these optimizers should be bound: gqcpy
 */
void bindDOCINewtonOrbitalOptimizers(py::module& module) {

    bindDOCINewtonOrbitalOptimizer<Algorithm<EigenproblemEnvironment>>(module, "Dense");
    bindDOCINewtonOrbitalOptimizer<IterativeAlgorithm<EigenproblemEnvironment>>(module, "Iterative");
}


}  // namespace gqcpy
