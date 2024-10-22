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

#include "Basis/NonOrthogonalBasis/GNonOrthogonalStateBasis.hpp"
#include "gqcpy/include/interfaces.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>


namespace py = pybind11;


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


/**
 *  Add Python bindings for some APIs related to `GNonOrthogonalStateBasis`.
 *
 *  @tparam Class               The type of the Pybind11 `class_` (generated by the compiler).
 *
 *  @param py_class             The Pybind11 `class_` that should obtain APIs related to `GNonOrthogonalStateBasis`.
 */
template <typename Class>
void bindGNonOrthogonalStateBasisInterface(Class& py_class) {

    // The C++ type corresponding to the Python class.
    using Type = typename Class::type;
    using Scalar = typename Type::Scalar;


    /**
     *  MARK: Constructors
     */

    py_class
        .def(py::init<const std::vector<GTransformation<Scalar>>&, const ScalarGSQOneElectronOperator<Scalar>&, const size_t, const double>(),
             py::arg("basis_state_vector"),
             py::arg("overlap_operator"),
             py::arg("number_of_occupied_orbitals"),
             py::arg("zero_threshold") = 1e-8);


    // Expose some Mulliken API to the Python class;
    bindNonOrthogonalStateBasisInterface(py_class);
}


/**
 *  Register `GNonOrthogonalStateBasis_d` and `GNonOrthogonalStateBasis_cd` to the gqcpy module and expose parts of their C++ interfaces to Python.
 *
 *  @param module           The Pybind11 module in which the classes should be registered.
 */
void bindGNonOrthogonalStateBases(py::module& module) {

    // Define the Python class for `GNonOrthogonalStateBasis_d`.
    py::class_<GNonOrthogonalStateBasis<double>> py_GNonOrthogonalStateBasis_d {module, "GNonOrthogonalStateBasis_d", "A class that represents a real, non-orthogonal state basis, created from `generalized` states."};

    bindGNonOrthogonalStateBasisInterface(py_GNonOrthogonalStateBasis_d);


    // Define the Python class for `GNonOrthogonalStateBasis_cd`.
    py::class_<GNonOrthogonalStateBasis<complex>> py_GNonOrthogonalStateBasis_cd {module, "GNonOrthogonalStateBasis_cd", "A class that represents a complex, non-orthogonal state basis, created from `generalized` states."};

    bindGNonOrthogonalStateBasisInterface(py_GNonOrthogonalStateBasis_cd);
}


}  // namespace gqcpy
