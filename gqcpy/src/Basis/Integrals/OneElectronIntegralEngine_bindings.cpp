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

#include "Basis/Integrals/OneElectronIntegralEngine.hpp"
#include "Basis/Integrals/Primitive/FunctionalOneElectronPrimitiveIntegralEngine.hpp"
#include "Utilities/aliases.hpp"

#include <pybind11/pybind11.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


/**
 *  Register `FunctionalOneElectronIntegralEngine_d` and `FunctionalOneElectronIntegralEngine_cd` to the gqcpy module and expose a part of their C++ interface to Python.
 * 
 *  @param module           The Pybind11 module in which `FunctionalOneElectronIntegralEngine_d` and `FunctionalOneElectronIntegralEngine_cd` should be registered.
 */
void bindFunctionalOneElectronIntegralEngine(py::module& module) {

    // Define the Python class for `FunctionalOneElectronIntegralEngine_d`.
    py::class_<OneElectronIntegralEngine<FunctionalOneElectronPrimitiveIntegralEngine<double>>> py_FunctionalOneElectronIntegralEngine_d {module, "FunctionalOneElectronIntegralEngine_d", "A custom one-electron integral engine that can calculate real-valued one-electron integrals over Cartesian GTO shells according to a custom implementation."};

    py_FunctionalOneElectronIntegralEngine_d
        .def(py::init<const FunctionalOneElectronPrimitiveIntegralEngine<double>&>(),
             py::arg("primitive_engine"));


    // Define the Python class for `FunctionalOneElectronIntegralEngine_cd`.
    py::class_<OneElectronIntegralEngine<FunctionalOneElectronPrimitiveIntegralEngine<complex>>> py_FunctionalOneElectronIntegralEngine_cd {module, "FunctionalOneElectronIntegralEngine_cd", "A custom one-electron integral engine that can calculate complex-valued one-electron integrals over Cartesian GTO shells according to a custom implementation."};

    py_FunctionalOneElectronIntegralEngine_cd
        .def(py::init<const FunctionalOneElectronPrimitiveIntegralEngine<complex>&>(),
             py::arg("primitive_engine"));
}


}  // namespace gqcpy
