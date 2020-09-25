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

#include "Processing/DensityMatrices/SpinResolvedTwoDM.hpp"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>


namespace py = pybind11;


namespace gqcpy {


void bindSpinResolvedTwoDM(py::module& module) {
    py::class_<GQCP::SpinResolvedTwoDM<Scalar>>(module, "SpinResolvedTwoDM", "A class that represents a spin resolved two DM.")

        // PUBLIC METHODS

        .def(
            "alphaAlpha",
            &GQCP::SpinResolvedTwoDM<Scalar>::alphaAlpha,
            "Return the pure alpha component of the spin resolved two DM")

        .def(
            "alphaBeta",
            &GQCP::SpinResolvedTwoDM<Scalar>::alphaBeta,
            "Return the mixed alpha-beta component of the spin resolved two DM")

        .def(
            "betaAlpha",
            &GQCP::SpinResolvedTwoDM<Scalar>::betaAlpha,
            "Return the mixed beta-alpha component of the spin resolved two DM")

        .def(
            "betaBeta",
            &GQCP::SpinResolvedTwoDM<Scalar>::betaBeta,
            "Return the pure beta component of the spin resolved two DM")

        .def(
            "numberOfOrbitals",
            [](const SpinResolvedTwoDM<Scalar>& d, const GQCP::Spin sigma, const GQCP::Spin tau) {
                return d.numberOfOrbitals(sigma, tau);
            },
            "Return the number of orbitals (spinors or spin-orbitals, depending on the context) that are related to the sigma-tau part of the spin-resolved 2-DM.")

        .def(
            "spinSummed",
            [](const SpinResolvedTwoDM<Scalar>& d) {
                return d.spinSummed();
            },
            "Return the spin-summed (total) 2-DM, i.e. the sum of four spin parts.")

    // PUBLIC METHODS
}