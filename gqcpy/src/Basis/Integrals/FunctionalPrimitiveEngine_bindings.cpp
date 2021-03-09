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

#include "Basis/Integrals/FunctionalPrimitiveEngine.hpp"
#include "Utilities/aliases.hpp"
#include "gqcpy/include/interfaces.hpp"

#include <pybind11/functional.h>
#include <pybind11/pybind11.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


/**
 *  Register `FunctionalPrimitiveEngine_d` and `FunctionalPrimitiveEngine_cd` to the gqcpy module and expose a part of their C++ interface to Python.
 * 
 *  @param module           The Pybind11 module in which `FunctionalPrimitiveEngine_d` and `FunctionalPrimitiveEngine_cd` should be registered.
 */
void bindFunctionalPrimitiveEngine(py::module& module) {

    // Define the Python class for `FunctionalPrimitiveEngine_d`.
    py::class_<FunctionalPrimitiveEngine<double>> py_FunctionalPrimitiveEngine_d {module, "FunctionalPrimitiveEngine_d", "A custom primitive engine that can real-valued calculate one-electron integrals over Cartesian GTOs according to a custom implementation."};

    py_FunctionalPrimitiveEngine_d
        .def(py::init<const std::function<double(const CartesianGTO&, const CartesianGTO&)>&>(),
             py::arg("function"));


    // Define the Python class for `FunctionalPrimitiveEngine_cd`.
    py::class_<FunctionalPrimitiveEngine<complex>>
        py_FunctionalPrimitiveEngine_cd {module, "FunctionalPrimitiveEngine_cd", "A custom primitive engine that can calculate complex-valued one-electron integrals over Cartesian GTOs according to a custom implementation.."};

    py_FunctionalPrimitiveEngine_cd
        .def(py::init<const std::function<complex(const CartesianGTO&, const CartesianGTO&)>&>(),
             py::arg("function"));
}


}  // namespace gqcpy
