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

#include "ONVBasis/SpinUnresolvedOperatorString.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


void bindSpinUnresolvedOperatorString(py::module& module) {
    py::class_<SpinUnresolvedOperatorString>(module, "SpinUnresolvedOperatorString", "A spin-unresolved operator string.")

        // CONSTRUCTORS

        .def(py::init<const std::vector<size_t>&>(),
             py::arg("index_vector"))

        .def(py::init<const std::vector<size_t>&, const int>(),
             py::arg("index_vector"),
             py::arg("phase_factor"))

        // NAMED CONSTRUCTORS

        .def_static(
            "FromONV",
            &SpinUnresolvedOperatorString::FromONV,
            py::arg("onv"))

        // PUBLIC METHODS

        .def(
            "operatorIndices",
            &SpinUnresolvedOperatorString::operatorIndices,
            "Retrieve the operator indices from the `SpinUnresolvedOperatorString`.")

        .def(
            "phaseFactor",
            &SpinUnresolvedOperatorString::phaseFactor,
            "Retrieve the phase factor corresponding to the `SpinUnresolvedOperatorString`.")

        .def(
            "isZero",
            &SpinUnresolvedOperatorString::isZero,
            "Check whether the operator string in question will result in zero when applied to the wave function.")

        .def(
            "sort",
            &SpinUnresolvedOperatorString::sort,
            "Sort the operator string in ascending order and adjust its phase factor.")

        .def(
            "schmidtDecomposition",
            &SpinUnresolvedOperatorString::schmidtDecomposition,
            "Decomposition of the `SpinUnresolvedOperatorString` into two new operator strings: a system and an environment.",
            py::arg("partition"));
}


}  // namespace gqcpy
