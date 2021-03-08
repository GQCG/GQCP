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

#include "Basis/Integrals/McMurchieDavidsonCoefficient.hpp"

#include <pybind11/pybind11.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


/**
 *  Register `McMurchieDavidsonCoefficient` to the gqcpy module and expose a part of its C++ interface to Python.
 * 
 *  @param module           The Pybind11 module in which `McMurchieDavidsonCoefficient` should be registered.
 */
void bindMcMurchieDavidsonCoefficient(py::module& module) {

    // Define the Python class for `GTransformation_d`.
    py::class_<McMurchieDavidsonCoefficient<double>> py_McMurchieDavidsonCoefficient {module, "McMurchieDavidsonCoefficient", "An implementation of the McMurchie-Davidson expansion coefficients through recurrence relations."};

    py_McMurchieDavidsonCoefficient
        .def(py::init<const double, const double, const double, const double>(),
             py::arg("K"),
             py::arg("a"),
             py::arg("L"),
             py::arg("b"))

        .def("__call__",
             py::arg("i"),
             py::arg("j"),
             py::arg("t"),
             "Return the value for the McMurchie-Davidson expansion coefficient E^{i,j}_t.");
}


}  // namespace gqcpy
