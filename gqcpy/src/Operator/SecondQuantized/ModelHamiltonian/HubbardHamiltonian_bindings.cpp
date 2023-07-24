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

#include "Operator/SecondQuantized/ModelHamiltonian/HubbardHamiltonian.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


void bindHubbardHamiltonian(py::module& module) {
    py::class_<HubbardHamiltonian<double>>(module, "HubbardHamiltonian", "The Hubbard model Hamiltonian.")

        /*
         *  MARK: Constructors
         */

        .def(py::init<const HoppingMatrix<double>&, const double&, const double&>(),
             py::arg("H"),
             py::arg("U"),
             py::arg("mu") = 0.0)

        .def(py::init<const HoppingMatrix<double>&, std::vector<double>&, std::vector<double>&>(),
             py::arg("H"),
             py::arg("U_vector"),
             py::arg("mu_vector"))


        /*
         *  MARK: Integral access
         */

        .def(
            "core",
            [](const HubbardHamiltonian<double>& hamiltonian) {
                return hamiltonian.core();
            },
            "Return the core Hamiltonian (i.e. resulting from the Hubbard hopping operator) as a one-electron operator.")

        .def(
            "twoElectron",
            [](const HubbardHamiltonian<double>& hamiltonian) {
                return hamiltonian.twoElectron();
            },
            "Return the two-electron part of the Hamiltonian (resulting from the on-site repulsion) as a two-electron operator.")

        .def(
            "hoppingMatrix",
            [](const HubbardHamiltonian<double>& hamiltonian) {
                return hamiltonian.hoppingMatrix();
            },
            "Return the Hubbard Hamiltonian hopping matrix for this Hubbard model Hamiltonian.")

        .def(
            "oneElectronContributions",
            [](const HubbardHamiltonian<double>& hamiltonian) {
                return hamiltonian.oneElectronContributions();
            },
            "Return the one electron contributions as a matrix for this Hubbard model Hamiltonian.")

        .def(
            "onSitePotentialMatrix",
            [](const HubbardHamiltonian<double>& hamiltonian) {
                return hamiltonian.onSitePotentialMatrix();
            },
            "Return the on-site potentials as a matrix for this Hubbard model Hamiltonian.")

        .def(
            "onSiteRepulsionMatrix",
            [](const HubbardHamiltonian<double>& hamiltonian) {
                return hamiltonian.onSiteRepulsionMatrix();
            },
            "Return the on-site repulsion contributions as a matrix for this Hubbard model Hamiltonian.");
}


}  // namespace gqcpy
