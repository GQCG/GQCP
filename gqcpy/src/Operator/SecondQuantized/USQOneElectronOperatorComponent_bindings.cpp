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

#include "Operator/SecondQuantized/USQOneElectronOperatorComponent.hpp"
#include "Utilities/aliases.hpp"
#include "gqcpy/include/interfaces.hpp"

#include <pybind11/pybind11.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


/**
 *  Register multiple variants of `USQOneElectronOperatorComponent` to the gqcpy module and expose a part of its C++ interface to Python.
 * 
 *  @param module           The Pybind11 module in which the `USQOneElectronOperatorComponent`s should be registered.
 */
void bindUSQOneElectronOperatorComponent(py::module& module) {

    // Define Python classes related to `USQOneElectronOperatorComponent` and expose their interfaces.
    py::class_<ScalarUSQOneElectronOperatorComponent<double>> py_ScalarUSQOneElectronOperatorComponent_d {module, "ScalarUSQOneElectronOperatorComponent_d", "One of the spin components of a (real) unrestricted one-electron operator, i.e. either the alpha or beta part."};

    bindSimpleSQOneElectronOperatorInterface(py_ScalarUSQOneElectronOperatorComponent_d);
    bindScalarSQOneElectronOperatorParameterInterface(py_ScalarUSQOneElectronOperatorComponent_d);


    py::class_<ScalarUSQOneElectronOperatorComponent<complex>> py_ScalarUSQOneElectronOperatorComponent_cd {module, "ScalarUSQOneElectronOperatorComponent_cd", "One of the spin components of a (complex) unrestricted one-electron operator, i.e. either the alpha or beta part."};

    bindSimpleSQOneElectronOperatorInterface(py_ScalarUSQOneElectronOperatorComponent_cd);
    bindScalarSQOneElectronOperatorParameterInterface(py_ScalarUSQOneElectronOperatorComponent_cd);


    py::class_<VectorUSQOneElectronOperatorComponent<double>> py_VectorUSQOneElectronOperatorComponent_d {module, "VectorUSQOneElectronOperatorComponent_d", "One of the spin components of a (real) unrestricted one-electron operator, i.e. either the alpha or beta part."};

    bindSimpleSQOneElectronOperatorInterface(py_VectorUSQOneElectronOperatorComponent_d);


    py::class_<VectorUSQOneElectronOperatorComponent<complex>> py_VectorUSQOneElectronOperatorComponent_cd {module, "VectorUSQOneElectronOperatorComponent_cd", "One of the spin components of a (complex) unrestricted one-electron operator, i.e. either the alpha or beta part."};

    bindSimpleSQOneElectronOperatorInterface(py_VectorUSQOneElectronOperatorComponent_cd);
}


}  // namespace gqcpy
