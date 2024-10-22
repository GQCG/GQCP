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

#include "Operator/SecondQuantized/GSQOneElectronOperator.hpp"
#include "Operator/SecondQuantized/GSQTwoElectronOperator.hpp"
#include "Utilities/complex.hpp"
#include "gqcpy/include/interfaces.hpp"

#include <pybind11/operators.h>
#include <pybind11/pybind11.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


/**
 *  Register multiple variants of `GSQOneElectronOperator` to the gqcpy module and expose a part of its C++ interface to Python.
 *
 *  @param module           The Pybind11 module in which the `GSQOneElectronOperator`s should be registered.
 */
void bindGSQOneElectronOperator(py::module& module) {

    // Scalar one-electron operators (real).
    py::class_<ScalarGSQOneElectronOperator<double>> py_ScalarGSQOneElectronOperator_d {module, "ScalarGSQOneElectronOperator_d", "A (real) general(ized) one-electron operator, which is suited for expressing spin-dependent one-electron operators."};

    bindSimpleSQOneElectronOperatorInterface(py_ScalarGSQOneElectronOperator_d);
    bindScalarSQOneElectronOperatorParameterInterface(py_ScalarGSQOneElectronOperator_d);

    py_ScalarGSQOneElectronOperator_d
        .def(py::self * py::self);


    // Scalar one-electron operators (complex).
    py::class_<ScalarGSQOneElectronOperator<complex>> py_ScalarGSQOneElectronOperator_cd {module, "ScalarGSQOneElectronOperator_cd", "A (complex) general(ized) one-electron operator, which is suited for expressing spin-dependent one-electron operators."};

    bindSimpleSQOneElectronOperatorInterface(py_ScalarGSQOneElectronOperator_cd);
    bindScalarSQOneElectronOperatorParameterInterface(py_ScalarGSQOneElectronOperator_cd);

    py_ScalarGSQOneElectronOperator_cd
        .def(py::self * py::self);


    // Vector one-electron operators (real).
    py::class_<VectorGSQOneElectronOperator<double>> py_VectorGSQOneElectronOperator_d {module, "VectorGSQOneElectronOperator_d", "A (real) general(ized) one-electron operator, which is suited for expressing spin-dependent one-electron operators."};

    bindSimpleSQOneElectronOperatorInterface(py_VectorGSQOneElectronOperator_d);


    // Vector one-electron operators (complex).
    py::class_<VectorGSQOneElectronOperator<complex>> py_VectorGSQOneElectronOperator_cd {module, "VectorGSQOneElectronOperator_cd", "A (complex) general(ized) one-electron operator, which is suited for expressing spin-dependent one-electron operators."};

    bindSimpleSQOneElectronOperatorInterface(py_VectorGSQOneElectronOperator_cd);
}


}  // namespace gqcpy
