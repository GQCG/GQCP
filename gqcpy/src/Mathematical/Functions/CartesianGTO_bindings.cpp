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

#include "Mathematical/Functions/CartesianGTO.hpp"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


/**
 *  Register `CartesianGTO` to the gqcpy module and expose a part of its C++ interface to Python.
 * 
 *  @param module           The Pybind11 module in which `CartesianGTO` should be registered.
 */
void bindCartesianGTO(py::module& module) {

    // Define the Python class for `CartesianGTO`.
    py::class_<CartesianGTO> py_CartesianGTO {module, "CartesianGTO", "A Cartesian Gaussian-type orbital (GTO), often referred to as a Cartesian 'primitive'."};


    py_CartesianGTO
        .def("gaussianExponent",
             &CartesianGTO::gaussianExponent,
             "Return the Gaussian exponent for this Cartesian GTO.")

        .def("center",
             &CartesianGTO::center,
             "Return the center of this Cartesian GTO.")

        .def("cartesianExponents",
             &CartesianGTO::cartesianExponents,
             "Return the Cartesian exponents for this Cartesian GTO.");
}


}  // namespace gqcpy
