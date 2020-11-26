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

#include "DensityMatrix/SpinDensity1DM.hpp"
#include "gqcpy/include/interfaces.hpp"

#include <pybind11/pybind11.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


/**
 *  Register `SpinDensity1DM_d` to the gqcpy module and expose a part of its C++ interface to Python.
 * 
 *  @param module           The Pybind11 module in which `SpinDensity1DM_d` should be registered.
 */
void bindSpinDensity1DM(py::module& module) {

    // Define the Python class for `SpinDensity1DM`.
    py::class_<SpinDensity1DM<double>> py_SpinDensity1DM_d {module, "SpinDensity1DM_d", "A type used to represent a one-electron spin-density density matrix, i.e. the alpha density matrix minus the beta density matrix."};


    // Expose the `BasisTransformable` API to the Python class.
    bindBasisTransformableInterface(py_SpinDensity1DM_d);
}


}  // namespace gqcpy
