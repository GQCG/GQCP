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

#include "QCMethod/Geminals/PSEnvironment.hpp"

#include <pybind11/pybind11.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


void bindPSEnvironment(py::module& module) {

    auto module_ps_environment = module.def_submodule("PSEnvironment");

    module_ps_environment.def(
        "AP1roG",
        [](const RSQHamiltonian<double>& sq_hamiltonian, const AP1roGGeminalCoefficients& G_initial) {
            return PSEnvironment::AP1roG(sq_hamiltonian, G_initial);
        },
        py::arg("sq_hamiltonian"),
        py::arg("G_initial"),
        "Return an environment suitable for AP1roG calculations.");

    module_ps_environment.def(
        "AP1roG",
        [](const RSQHamiltonian<double>& sq_hamiltonian, const size_t N_P) {
            return PSEnvironment::AP1roG(sq_hamiltonian, N_P);
        },
        py::arg("sq_hamiltonian"),
        py::arg("N_P"),
        "Return an environment suitable for AP1roG calculations.");
}


}  // namespace gqcpy
