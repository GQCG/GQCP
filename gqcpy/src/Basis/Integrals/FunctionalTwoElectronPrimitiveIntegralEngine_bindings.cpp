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

#include "Basis/Integrals/Primitive/FunctionalTwoElectronPrimitiveIntegralEngine.hpp"
#include "Utilities/aliases.hpp"
#include "gqcpy/include/interfaces.hpp"

#include <pybind11/functional.h>
#include <pybind11/pybind11.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


/**
 *  Register `FunctionalTwoElectronPrimitiveIntegralEngine_d` and `FunctionalTwoElectronPrimitiveIntegralEngine_cd` to the gqcpy module and expose a part of their C++ interface to Python.
 * 
 *  @param module           The Pybind11 module in which `FunctionalTwoElectronPrimitiveIntegralEngine_d` and `FunctionalTwoElectronPrimitiveIntegralEngine_cd` should be registered.
 */
void bindFunctionalTwoElectronPrimitiveIntegralEngine(py::module& module) {

    // Define the Python class for `FunctionalTwoElectronPrimitiveIntegralEngine_d`.
    py::class_<FunctionalTwoElectronPrimitiveIntegralEngine<double>> py_FunctionalTwoElectronPrimitiveIntegralEngine_d {module, "FunctionalTwoElectronPrimitiveIntegralEngine_d", "A custom primitive engine that can calculate real-valued two-electron integrals over Cartesian GTOs according to a custom implementation."};

    py_FunctionalTwoElectronPrimitiveIntegralEngine_d
        .def(py::init<const std::function<double(const CartesianGTO&, const CartesianGTO&, const CartesianGTO&, const CartesianGTO&)>&>(),
             py::arg("function"));


    // Define the Python class for `FunctionalTwoElectronPrimitiveIntegralEngine_cd`.
    py::class_<FunctionalTwoElectronPrimitiveIntegralEngine<complex>>
        py_FunctionalTwoElectronPrimitiveIntegralEngine_cd {module, "FunctionalTwoElectronPrimitiveIntegralEngine_cd", "A custom primitive engine that can calculate complex-valued two-electron integrals over Cartesian GTOs according to a custom implementation."};

    py_FunctionalTwoElectronPrimitiveIntegralEngine_cd
        .def(py::init<const std::function<complex(const CartesianGTO&, const CartesianGTO&, const CartesianGTO&, const CartesianGTO&)>&>(),
             py::arg("function"));
}


}  // namespace gqcpy
