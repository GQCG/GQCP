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

#include "Basis/ScalarBasis/GTOShell.hpp"
#include "Basis/SpinorBasis/RSpinOrbitalBasis.hpp"
#include "Molecule/Molecule.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>


namespace py = pybind11;


namespace gqcpy {


void bindSQHamiltonian(py::module& module) {
    py::class_<GQCP::RSQHamiltonian<double>>(module, "RSQHamiltonian", "A class that represents a real, restricted second-quantized Hamiltonian.")

        // CONSTRUCTORS

        .def_static(
            "Molecular",
            [](const GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell>& r_spinor_basis, const GQCP::Molecule& molecule) {
                return GQCP::RSQHamiltonian<double>::Molecular(r_spinor_basis, molecule);
            },
            py::arg("r_spinor_basis"),
            py::arg("molecule"),
            "Construct the molecular Hamiltonian in a given restricted spin-orbital basis.")

        // .def_static(
        //     "Molecular",
        //     [](const GQCP::GSpinorBasis<double, GQCP::GTOShell>& g_spinor_basis, const GQCP::Molecule& molecule) {
        //         return GQCP::RSQHamiltonian<double>::Molecular(g_spinor_basis, molecule);
        //     },
        //     py::arg("g_spinor_basis"),
        //     py::arg("molecule"),
        //     "Construct the molecular Hamiltonian in a given (general) spinor basis.")


        // PUBLIC METHODS

        .def(
            "__add__",
            [](const GQCP::RSQHamiltonian<double>& sq_hamiltonian, const GQCP::ScalarRSQOneElectronOperator<double>& sq_op) {
                return sq_hamiltonian + sq_op;
            })

        .def(
            "__sub__",
            [](const GQCP::RSQHamiltonian<double>& sq_hamiltonian, const GQCP::ScalarRSQOneElectronOperator<double>& sq_op) {
                return sq_hamiltonian - sq_op;
            });
}


}  // namespace gqcpy
