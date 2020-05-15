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

#include "QCMethod/Geminals/AP1roGLagrangianNewtonOrbitalOptimizer.hpp"

#include "Mathematical/Optimization/Minimization/IterativeIdentitiesHessianModifier.hpp"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>


namespace py = pybind11;


namespace gqcpy {


void bindAP1roGLagrangianNewtonOrbitalOptimizer(py::module& module) {
    py::class_<GQCP::AP1roGLagrangianNewtonOrbitalOptimizer>(module,
                                                             "AP1roGLagrangianNewtonOrbitalOptimizer",
                                                             "An orbital optimizer for vAP1roG.")

        // Use a standard Hessian modifier for the Python bindings.
        .def(py::init([](const GQCP::AP1roGGeminalCoefficients& G, const double oo_convergence_threshold = 1.0e-08, const size_t oo_maximum_number_of_iterations = 128, const double pse_convergence_threshold = 1.0e-08, const size_t pse_maximum_number_of_iterations = 128) {
                 auto hessian_modifier = std::make_shared<GQCP::IterativeIdentitiesHessianModifier>();
                 return GQCP::AP1roGLagrangianNewtonOrbitalOptimizer(G, hessian_modifier, oo_convergence_threshold, oo_maximum_number_of_iterations, pse_convergence_threshold, pse_maximum_number_of_iterations);
             }),
             py::arg("G"),
             py::arg("oo_convergence_threshold") = 1.0e-08,
             py::arg("oo_maximum_number_of_iterations") = 128,
             py::arg("pse_convergence_threshold") = 1.0e-08,
             py::arg("pse_maximum_number_of_iterations") = 128)

        .def("get_electronic_energy",
             &GQCP::AP1roGLagrangianNewtonOrbitalOptimizer::get_electronic_energy)

        .def("get_geminal_coefficients",
             &GQCP::AP1roGLagrangianNewtonOrbitalOptimizer::get_geminal_coefficients)

        .def("get_multipliers",
             [](const GQCP::AP1roGLagrangianNewtonOrbitalOptimizer& optimizer) {
                 return optimizer.get_multipliers().asMatrix();
             })

        .def("optimize",
             [](GQCP::AP1roGLagrangianNewtonOrbitalOptimizer& optimizer, GQCP::RSpinorBasis<double, GQCP::GTOShell>& spinor_basis, GQCP::SQHamiltonian<double>& sq_hamiltonian) {
                 return optimizer.optimize(spinor_basis, sq_hamiltonian);
             });
}


}  // namespace gqcpy
