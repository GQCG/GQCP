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

#include "Basis/Transformations/UTransformation.hpp"
#include "gqcpy/include/interfaces.hpp"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


/**
 *  Register `UTransformation_d` to the gqcpy module and expose a part of its C++ interface to Python.
 * 
 *  @param module           The Pybind11 module in which `UTransformation_d` should be registered.
 */
void bindUTransformation(py::module& module) {

    // Define the Python class for `UTransformation`.
    py::class_<UTransformation<double>> py_UTransformation_d {module, "UTransformation_d", "A type that encapsulates transformation matrices for the alpha- and beta-parts of spin-orbital bases."};

    py_UTransformation_d
        .def("inverse",
             &UTransformation<double>::inverse,
             "Return the inverse transformation of this transformation matrix.")

        .def(
            "isUnitary",
            [](const UTransformation<double>& T, const double threshold = 1.0e-12) {
                return T.isUnitary(threshold);
            },
            "Return if this transformation is unitary, within the given threshold");

    // Expose the `SpinResolvedBase` API to Python.
    bindSpinResolvedBaseInterface(py_UTransformation_d);

    // Expose the `BasisTransformable` API to Python.
    bindBasisTransformableInterface(py_UTransformation_d);
}


}  // namespace gqcpy
