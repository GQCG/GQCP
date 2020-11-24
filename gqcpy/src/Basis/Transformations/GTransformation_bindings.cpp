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

#include "Basis/Transformations/GTransformation.hpp"
#include "gqcpy/include/interfaces.hpp"

#include <pybind11/pybind11.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


/**
 *  Register `GTransformation_d` to the gqcpy module and expose a part of its C++ interface to Python.
 * 
 *  @param module           The Pybind11 module in which `UTransformation_d` should be registered.
 */
void bindGTransformation(py::module& module) {

    // Define the Python class for `GTransformation_d`.
    py::class_<GTransformation<double>> py_GTransformation_d {module, "GTransformation_d", "A 'general' basis transformation, i.e. a general, full-spinor basis transformation where the transformation mixes the alpha- and beta components of the two-component spinors."};


    // Expose the `SimpleTransformation` API to the Python class.
    bindSimpleTransformationInterface(py_GTransformation_d);
}


}  // namespace gqcpy
