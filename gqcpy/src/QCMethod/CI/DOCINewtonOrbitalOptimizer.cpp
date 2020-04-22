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
#include "QCMethod/CI/DOCINewtonOrbitalOptimizer.hpp"

#include "Mathematical/Algorithm/Algorithm.hpp"
#include "Mathematical/Algorithm/IterativeAlgorithm.hpp"

#include <pybind11/pybind11.h>


namespace py = pybind11;


namespace gqcpy {


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

    py::class_<GQCP::DOCINewtonOrbitalOptimizer<EigenproblemSolver>>(module,
                                                                     ("DOCINewtonOrbitalOptimizer_" + suffix).c_str(),
                                                                     "A class that performs gradient-and-Hessian-based orbital optimization for DOCI")

        .def(
            "optimize",
            [](GQCP::DOCINewtonOrbitalOptimizer<EigenproblemSolver>& optimizer, GQCP::RSpinorBasis<double, GQCP::GTOShell>& spinor_basis, GQCP::SQHamiltonian<double>& sq_hamiltonian) {
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

    bindDOCINewtonOrbitalOptimizer<GQCP::Algorithm<GQCP::EigenproblemEnvironment>>(module, "Dense");
    bindDOCINewtonOrbitalOptimizer<GQCP::IterativeAlgorithm<GQCP::EigenproblemEnvironment>>(module, "Iterative");
}


}  // namespace gqcpy
