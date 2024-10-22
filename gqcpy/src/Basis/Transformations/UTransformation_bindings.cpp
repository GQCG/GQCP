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

    // Define the Python class for `UTransformation_d`.
    py::class_<UTransformation<double>> py_UTransformation_d {module, "UTransformation_d", "A type that encapsulates transformation matrices for the alpha- and beta-parts of spin-orbital bases."};

    py_UTransformation_d

        /*
         *  MARK: Named constructors
         */

        .def_static(
            "FromRestricted",
            [](const RTransformation<double>& T) {
                return UTransformation<double>::FromRestricted(T);
            },
            py::arg("T"),
            "Create an unrestricted transformation from an restricted transformation.");


    // Expose the `SpinResolvedBase` API to Python.
    bindSpinResolvedBaseInterface(py_UTransformation_d);

    // Expose the `BasisTransformable` API to Python.
    bindBasisTransformableInterface(py_UTransformation_d);

    // Add some APIs related to operations on `BasisTransformable` objects.
    bindBasisTransformableOperationsInterface(py_UTransformation_d);

    // Define the Python class for `UTransformation_cd`.
    py::class_<UTransformation<complex>> py_UTransformation_cd {module, "UTransformation_cd", "A type that encapsulates transformation matrices for the alpha- and beta-parts of spin-orbital bases."};


    py_UTransformation_cd

        /*
         *  MARK: Named constructors
         */

        .def_static(
            "FromRestricted",
            [](const RTransformation<complex>& T) {
                return UTransformation<complex>::FromRestricted(T);
            },
            py::arg("T"),
            "Create an unrestricted transformation from an restricted transformation.");

    // Expose the `SpinResolvedBase` API to Python.
    bindSpinResolvedBaseInterface(py_UTransformation_cd);

    // Expose the `BasisTransformable` API to Python.
    bindBasisTransformableInterface(py_UTransformation_cd);

    // Add some APIs related to operations on `BasisTransformable` objects.
    bindBasisTransformableOperationsInterface(py_UTransformation_cd);
}


}  // namespace gqcpy
