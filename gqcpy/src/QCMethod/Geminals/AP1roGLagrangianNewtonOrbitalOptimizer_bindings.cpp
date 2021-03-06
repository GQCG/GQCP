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
#include "QCMethod/Geminals/AP1roGLagrangianNewtonOrbitalOptimizer.hpp"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


void bindAP1roGLagrangianNewtonOrbitalOptimizer(py::module& module) {
    py::class_<AP1roGLagrangianNewtonOrbitalOptimizer>(module,
                                                       "AP1roGLagrangianNewtonOrbitalOptimizer",
                                                       "An orbital optimizer for vAP1roG.")

        // CONSTRUCTORS

        // Use a standard Hessian modifier for the Python bindings.
        .def(py::init([](const AP1roGGeminalCoefficients& G, const double oo_convergence_threshold = 1.0e-08, const size_t oo_maximum_number_of_iterations = 128, const double pse_convergence_threshold = 1.0e-08, const size_t pse_maximum_number_of_iterations = 128) {
                 auto hessian_modifier = std::make_shared<IterativeIdentitiesHessianModifier>();
                 return AP1roGLagrangianNewtonOrbitalOptimizer(G, hessian_modifier, oo_convergence_threshold, oo_maximum_number_of_iterations, pse_convergence_threshold, pse_maximum_number_of_iterations);
             }),
             py::arg("G"),
             py::arg("oo_convergence_threshold") = 1.0e-08,
             py::arg("oo_maximum_number_of_iterations") = 128,
             py::arg("pse_convergence_threshold") = 1.0e-08,
             py::arg("pse_maximum_number_of_iterations") = 128)


        // PUBLIC METHODS

        .def("electronicEnergy",
             &AP1roGLagrangianNewtonOrbitalOptimizer::electronicEnergy)

        .def("geminalCoefficients",
             &AP1roGLagrangianNewtonOrbitalOptimizer::geminalCoefficients)

        .def("multipliers",
             [](const AP1roGLagrangianNewtonOrbitalOptimizer& optimizer) {
                 return optimizer.multipliers().asMatrix();
             })

        .def("optimize",
             [](AP1roGLagrangianNewtonOrbitalOptimizer& optimizer, RSpinOrbitalBasis<double, GTOShell>& spinor_basis, RSQHamiltonian<double>& sq_hamiltonian) {
                 return optimizer.optimize(spinor_basis, sq_hamiltonian);
             });
}


}  // namespace gqcpy
