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

#include "Processing/DensityMatrices/TwoDM.hpp"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>


namespace py = pybind11;


namespace gqcpy {


void bindTwoDM(py::module& module) {

    py::class_<GQCP::TwoDM<double>>(module, "TwoDM", "A two-particle density matrix")

        // PUBLIC METHODS

        .def(
            "reduce",
            [](const GQCP::TwoDM<double>& d) {
                return d.reduce();
            },
            "Return a partial contraction of the 2-DM, where D(p,q) = d(p,q,r,r).")

        .def(
            "trace",
            [](const GQCP::TwoDM<double>& d) {
                return d.trace();
            },
            "Return the trace of the 2-DM, i.e. d(p,p,q,q).");
}


}  // namespace gqcpy