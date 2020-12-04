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

#include "QuantumChemical/SpinResolved.hpp"
#include "gqcpy/include/interfaces.hpp"

#include <pybind11/pybind11.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;

template <typename Of>
void bindSpinResolved(py::module& module, const std::string& name, const std::string& description) {

    // Define the Python class for `spinResolved` types.
    py::class_<SpinResolved<Of>> py_spinResolved {module, name.c_str(), description.c_str()};

    // Expose the `BasisTransformable` API to the Python class.
    bindSpinResolvedBaseInterface(py_spinResolved);
}


/**
 *  Bind all types of `SpinResolved`s.
 */
void bindSpinResolvedTypes(py::module& module) {

    bindSpinResolved<size_t>(module, "SpinResolvedLong", "A spin resolved encapsulation of two unsigned longs.");
    bindSpinResolved<std::vector<double>>(module, "SpinResolvedStdVector", "A spin resolved encapsulation of two std::vectors.");
    bindSpinResolved<VectorX<double>>(module, "SpinResolvedVector", "A spin resolved encapsulation of two GQCP::Vectors.");
}

}  // namespace gqcpy
