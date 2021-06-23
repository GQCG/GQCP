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
#include "Basis/ScalarBasis/ScalarBasis.hpp"

#include <pybind11/pybind11.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


/**
 *  Register `ScalarBasis` to the gqcpy module and expose a part of its C++ interface to Python.
 * 
 *  @param module           The Pybind11 module in which `ScalarBasis` should be registered.
 */
void bindScalarBasis(py::module& module) {
    py::class_<ScalarBasis<GTOShell>>(module, "ScalarBasis", "A class that represents a scalar basis: it represents a collection of scalar basis functions.")

        .def("shellSet",
             &ScalarBasis<GTOShell>::shellSet,
             "The collection of shells that represents this scalar basis.");
}


}  // namespace gqcpy
