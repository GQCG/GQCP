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

#include "Mathematical/Optimization/Minimization/IterativeIdentitiesHessianModifier.hpp"
#include "QCMethod/CI/DOCINewtonOrbitalOptimizer.hpp"

#include <pybind11/pybind11.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


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
        [](const SeniorityZeroONVBasis& onv_basis, const EigenproblemSolver& solver, const EigenproblemEnvironment& environment, const size_t number_of_requested_eigenpairs = 1, const double convergence_threshold = 1.0e-08, const size_t maximum_number_of_iterations = 128) {
            auto hessian_modifier = std::make_shared<IterativeIdentitiesHessianModifier>();

            return DOCINewtonOrbitalOptimizer<EigenproblemSolver>(onv_basis, solver, environment, hessian_modifier, number_of_requested_eigenpairs, convergence_threshold, maximum_number_of_iterations);
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

    bindDOCINewtonOrbitalOptimizerFactoryMethod<IterativeAlgorithm<EigenproblemEnvironment>>(module);
    bindDOCINewtonOrbitalOptimizerFactoryMethod<Algorithm<EigenproblemEnvironment>>(module);
}


}  // namespace gqcpy
