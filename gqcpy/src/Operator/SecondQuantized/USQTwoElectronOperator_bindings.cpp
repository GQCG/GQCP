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

#include "Operator/SecondQuantized/USQTwoElectronOperator.hpp"
#include "gqcpy/include/interfaces.hpp"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>


namespace py = pybind11;


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


/**
 *  Register `ScalarUSQTwoElectronOperator_d` to the gqcpy module and expose a part of its C++ interface to Python.
 * 
 *  @param module           The Pybind11 module in which `ScalarUSQTwoElectronOperator_d` should be registered.
 */
void bindUSQTwoElectronOperator(py::module& module) {

    // Define the Python class for `USQTwoElectronOperator`.
    py::class_<ScalarUSQTwoElectronOperator<double>> py_ScalarUSQTwoElectronOperator_d {module, "ScalarUSQTwoElectronOperator_d", "A class that represents a (real) 'unrestricted second-quantized two-electron operator' suitable for the projection of the non-relativistic Hamiltonian onto an unrestricted spinor basis. It holds the tensor representation of its parameters for both spin components and both mixed spin components, which are (usually) integrals over first-quantized operators."};

    // Expose the `DoublySpinResolvedBase` API to the Python class.
    bindDoublySpinResolvedBaseInterface(py_ScalarUSQTwoElectronOperator_d);
}


}  // namespace gqcpy
