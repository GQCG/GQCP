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

#include "Operator/SecondQuantized/GSQTwoElectronOperator.hpp"
#include "Utilities/aliases.hpp"
#include "gqcpy/include/interfaces.hpp"

#include <pybind11/pybind11.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


/**
 *  Register multiple variants of `GSQTwoElectronOperator` to the gqcpy module and expose a part of its C++ interface to Python.
 * 
 *  @param module           The Pybind11 module in which the `GSQTwoElectronOperator`s should be registered.
 */
void bindGSQTwoElectronOperator(py::module& module) {

    // Define Python classes related to `GSQTwoElectronOperator` and expose their interfaces.
    py::class_<ScalarGSQTwoElectronOperator<double>> py_ScalarGSQTwoElectronOperator_d {module, "ScalarGSQTwoElectronOperator_d", "A (real) general(ized) two-electron operator, which is suited for expressing spin-dependent two-electron operators."};

    bindSimpleSQTwoElectronOperatorInterface(py_ScalarGSQTwoElectronOperator_d);
    bindScalarSQTwoElectronOperatorParameterInterface(py_ScalarGSQTwoElectronOperator_d);


    py::class_<ScalarGSQTwoElectronOperator<complex>> py_ScalarGSQTwoElectronOperator_cd {module, "ScalarGSQTwoElectronOperator_cd", "A (complex) general(ized) two-electron operator, which is suited for expressing spin-dependent two-electron operators."};

    bindSimpleSQTwoElectronOperatorInterface(py_ScalarGSQTwoElectronOperator_cd);
    bindScalarSQTwoElectronOperatorParameterInterface(py_ScalarGSQTwoElectronOperator_cd);
}


}  // namespace gqcpy
