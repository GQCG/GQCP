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

#include "Processing/Properties/DysonOrbital.hpp"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>


namespace py = pybind11;


namespace gqcpy {


void bindDysonOrbital(py::module& module) {

    py::class_<GQCP::DysonOrbital<double>>(module, "DysonOrbital", "A class representing a Dyson orbital expanded in a basis of spin-orbitals.")

        // NAMED CONSTRUCTORS

        .def_static(
            "FromOverlap",
            [](const GQCP::LinearExpansion<GQCP::SpinResolvedONVBasis>& linear_expansion1, const GQCP::LinearExpansion<GQCP::SpinResolvedONVBasis>& linear_expansion2) {
                return GQCP::DysonOrbital<double>::FromOverlap(linear_expansion1, linear_expansion2);
            },
            py::arg("linear_expansion1"),
            py::arg("linear_expansion2"),
            "Return a Dyson orbital containing Dyson amplitudes.")


        // PUBLIC METHODS

        .def(
            "amplitudes",
            &GQCP::DysonOrbital<double>::amplitudes,
            "Dyson amplitudes which are the expansion coefficients in a spin-orbital basis.");
}


}  // namespace gqcpy
