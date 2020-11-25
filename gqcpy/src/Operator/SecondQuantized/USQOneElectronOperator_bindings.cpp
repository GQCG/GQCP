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

#include "Basis/Transformations/USQOneElectronOperator.hpp"
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

    // Define the Python class for `USQOneElectronOperator`.
    py::class_<USQOneElectronOperator<double>> py_USQOneElectronOperator_d {module, "USQOneElectronOperator_d", "A class that represents a (real) 'unrestricted second-quantized one-electron operator'. This type of operator is suitable for the projection of the non-relativistic Hamiltonian onto an unrestricted spinor basis. It holds the matrix representation of its parameters for both spin components."};


    // Expose the `SpinResolvedBase` API to Python.
    bindSpinResolvedBaseInterface(py_USQOneElectronOperator_d);

    // Expose one-electron operator APIs to Python.
    bindSQOneElectronOperatorInterface(py_USQOneElectronOperator_d);
}


}  // namespace gqcpy
