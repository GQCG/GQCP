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

#include "ONVBasis/SpinResolvedOperatorString.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


void bindSpinResolvedOperatorString(py::module& module) {
    py::class_<SpinResolvedOperatorString>(module, "SpinResolvedOperatorString", "A spin-resolved operator string.")

        // CONSTRUCTORS

        .def(py::init<const std::vector<size_t>&, const std::vector<Spin>&>(),
             py::arg("index_vector"),
             py::arg("spin_vector"))

        // NAMED CONSTRUCTORS

        .def_static(
            "FromONV",
            &SpinResolvedOperatorString::FromONV,
            py::arg("onv"))

        // PUBLIC METHODS

        .def(
            "operatorIndices",
            &SpinResolvedOperatorString::operatorIndices,
            "Retrieve the operator indices from the `SpinResolvedOperatorString`.")

        .def(
            "operatorSpins",
            &SpinResolvedOperatorString::operatorSpins,
            "Retrieve the operator spins from the `SpinResolvedOperatorString`.")

        .def(
            "size",
            &SpinResolvedOperatorString::size,
            "Retrieve the number of operators in the `SpinResolvedOperatorString`.")

        .def(
            "spinResolve",
            &SpinResolvedOperatorString::spinResolve,
            "Return the spin resolved operator string split in its two separate components, one containing the alpha operators, one containing the beta operators. The alpha operator string will recieve the total phase factor assiciated with this action, the phase factor of the beta string will be one.");
}


}  // namespace gqcpy
