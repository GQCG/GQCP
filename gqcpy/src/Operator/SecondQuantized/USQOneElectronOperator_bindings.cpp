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

#include "Operator/SecondQuantized/USQOneElectronOperator.hpp"
#include "gqcpy/include/interfaces.hpp"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


/**
 *  Register `USQOneElectronOperator_d` to the gqcpy module and expose a part of its C++ interface to Python.
 * 
 *  @param module           The Pybind11 module in which `USQOneElectronOperator_d` should be registered.
 */
void bindUSQOneElectronOperator(py::module& module) {

    // Define the Python classes for `USQOneElectronOperator` and expose its APIs.
    py::class_<ScalarUSQOneElectronOperator<double>> py_ScalarUSQOneElectronOperator_d {module, "USQOneElectronOperator_d", "A class that represents a (real) 'unrestricted second-quantized one-electron operator'. This type of operator is suitable for the projection of the non-relativistic Hamiltonian onto an unrestricted spinor basis. It holds the matrix representation of its parameters for both spin components."};

    bindSpinResolvedBaseInterface(py_ScalarUSQOneElectronOperator_d);
    bindSQOneElectronOperatorInterface(py_ScalarUSQOneElectronOperator_d);


    py::class_<VectorUSQOneElectronOperator<double>> py_VectorUSQOneElectronOperator_d {module, "VectorUSQOneElectronOperator_d", "A class that represents a (real) 'unrestricted second-quantized one-electron operator'. This type of operator is suitable for the projection of the non-relativistic Hamiltonian onto an unrestricted spinor basis. It holds the matrix representation of its parameters for both spin components."};

    bindSpinResolvedBaseInterface(py_VectorUSQOneElectronOperator_d);
    bindSQOneElectronOperatorInterface(py_VectorUSQOneElectronOperator_d);
}


}  // namespace gqcpy
