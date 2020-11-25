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

#include "Operator/SecondQuantized/MixedUSQTwoElectronOperatorComponent.hpp"
#include "Utilities/aliases.hpp"
#include "gqcpy/include/interfaces.hpp"

#include <pybind11/pybind11.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


/**
 *  Register multiple variants of `MixedUSQTwoElectronOperatorComponent` to the gqcpy module and expose a part of its C++ interface to Python.
 * 
 *  @param module           The Pybind11 module in which the `MixedUSQTwoElectronOperatorComponent`s should be registered.
 */
void bindMixedUSQTwoElectronOperatorComponent(py::module& module) {

    // Define Python classes related to `MixedUSQTwoElectronOperatorComponent` and expose their interfaces.
    py::class_<ScalarMixedUSQTwoElectronOperatorComponent<double>> py_ScalarMixedUSQTwoElectronOperatorComponent_d {module, "ScalarMixedUSQTwoElectronOperatorComponent_d", "One of the mixed (i.e. alpha-beta or beta-alpha) spin components of a (real) unrestricted two-electron operator."};

    bindSQTwoElectronOperatorInterface(py_ScalarMixedUSQTwoElectronOperatorComponent_d);
    bindScalarSQTwoElectronOperatorParameterInterface(py_ScalarMixedUSQTwoElectronOperatorComponent_d);


    py::class_<ScalarMixedUSQTwoElectronOperatorComponent<complex>> py_ScalarMixedUSQTwoElectronOperatorComponent_cd {module, "ScalarMixedUSQTwoElectronOperatorComponent_cd", "One of the mixed (i.e. alpha-beta or beta-alpha) spin components of a (complex) unrestricted two-electron operator."};

    bindSQTwoElectronOperatorInterface(py_ScalarMixedUSQTwoElectronOperatorComponent_cd);
    bindScalarSQTwoElectronOperatorParameterInterface(py_ScalarMixedUSQTwoElectronOperatorComponent_cd);
}


}  // namespace gqcpy
