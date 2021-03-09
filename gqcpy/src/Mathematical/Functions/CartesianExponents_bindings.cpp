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

#include "Mathematical/Functions/CartesianExponents.hpp"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


/**
 *  Register `CartesianExponents` to the gqcpy module and expose a part of its C++ interface to Python.
 * 
 *  @param module           The Pybind11 module in which `CartesianExponents` should be registered.
 */
void bindCartesianExponents(py::module& module) {

    // Define the Python class for `CartesianExponents`.
    py::class_<CartesianExponents> py_CartesianExponents {module, "CartesianExponents", "A class that represents exponents of the Cartesian functions x, y and z."};

    py_CartesianExponents
        .def("value",
             &CartesianExponents::value,
             py::arg("direction"),
             "Return the exponent in the given direction");
}


}  // namespace gqcpy
