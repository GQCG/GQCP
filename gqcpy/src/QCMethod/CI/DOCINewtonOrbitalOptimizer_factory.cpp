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
#include "Mathematical/Optimization/Minimization/IterativeIdentitiesHessianModifier.hpp"
#include "QCMethod/CI/DOCINewtonOrbitalOptimizer.hpp"

#include <pybind11/pybind11.h>


namespace py = pybind11;


namespace gqcpy {


/**
 *  We want to provide an easier Python API for the construction of a DOCI newton-orbital optimizer methods. The following bindings implement this kind of factory behaviour.
 */

/**
 *  Bind a factory-like method for a DOCI Newton orbital optimizer.
 * 
 *  @tparam EigenproblemSolver          the type of the eigenvalue problem solver that should be used
 */
template <typename EigenproblemSolver>
void bindDOCINewtonOrbitalOptimizerFactoryMethod(py::module& module) {

    module.def(
        "DOCINewtonOrbitalOptimizer",
        [](const GQCP::SeniorityZeroONVBasis& onv_basis, const EigenproblemSolver& solver, const GQCP::EigenproblemEnvironment& environment, const size_t number_of_requested_eigenpairs = 1, const double convergence_threshold = 1.0e-08, const size_t maximum_number_of_iterations = 128) {
            auto hessian_modifier = std::make_shared<GQCP::IterativeIdentitiesHessianModifier>();

            return GQCP::DOCINewtonOrbitalOptimizer<EigenproblemSolver>(onv_basis, solver, environment, hessian_modifier, number_of_requested_eigenpairs, convergence_threshold, maximum_number_of_iterations);
        },
        py::arg("onv_basis"),
        py::arg("solver"),
        py::arg("environment"),
        py::arg("number_of_requested_eigenpairs") = 1,
        py::arg("convergence_threshold") = 1.0e-08,
        py::arg("maximum_number_of_iterations") = 128);
}


/**
 *  Bind all types of DOCI Newton orbital optimizers to the given module.
 */
void bindDOCINewtonOrbitalOptimizerFactory(py::module& module) {

    bindDOCINewtonOrbitalOptimizerFactoryMethod<GQCP::IterativeAlgorithm<GQCP::EigenproblemEnvironment>>(module);
    bindDOCINewtonOrbitalOptimizerFactoryMethod<GQCP::Algorithm<GQCP::EigenproblemEnvironment>>(module);
}


}  // namespace gqcpy
