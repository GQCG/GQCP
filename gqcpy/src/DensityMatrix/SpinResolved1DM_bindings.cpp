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

#include "DensityMatrix/SpinResolved1DM.hpp"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>


namespace py = pybind11;


namespace gqcpy {


void bindSpinResolvedOneDM(py::module& module) {
    py::class_<GQCP::SpinResolved1DM<double>>(module, "SpinResolved1DM", "A class that represents a spin resolved 1-DM.")

        // CONSTRUCTORS

        .def_static(
            "FromOrbital1DM",
            [](const Eigen::MatrixXd& D) {
                return GQCP::SpinResolved1DM<double>::FromOrbital1DM(GQCP::Orbital1DM<double> {D});
            },
            "Return a spin resolved 1-DM created from an orbital 1-DM.")

        // PUBLIC METHODS

        .def(
            "alpha",
            [](const GQCP::SpinResolved1DM<double>& D) {
                return D.alpha();
            },
            "Return the alpha part of the spin resolved 1-DM.")

        .def(
            "beta",
            [](const GQCP::SpinResolved1DM<double>& D) {
                return D.beta();
            },
            "Return the beta part of the spin resolved 1-DM.")

        .def(
            "numberOfOrbitals",
            [](const GQCP::SpinResolved1DM<double>& D, const GQCP::Spin sigma) {
                return D.numberOfOrbitals(sigma);
            },
            py::arg("sigma"),
            "Return the number of orbitals (spinors or spin-orbitals, depending on the context) that correspond to the given spin.")

        .def(
            "spinDensity",
            &GQCP::SpinResolved1DM<double>::spinDensity,
            "Return the spin-density matrix, i.e. the difference between the alpha and beta 1-DM.")

        .def(
            "orbitalDensity",
            &GQCP::SpinResolved1DM<double>::orbitalDensity,
            "Return the orbital density matrix, i.e. the sum of the alpha and beta 1-DM.");
}

}  // namespace gqcpy