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

#include "Basis/Integrals/HermiteCoulombIntegral.hpp"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


/**
 *  Register `HermiteCoulombIntegral` to the gqcpy module and expose a part of its C++ interface to Python.
 * 
 *  @param module           The Pybind11 module in which `HermiteCoulombIntegral` should be registered.
 */
void bindHermiteCoulombIntegral(py::module& module) {

    // Define the Python class for `HermiteCoulombIntegral`.
    py::class_<HermiteCoulombIntegral> py_HermiteCoulombIntegral {module, "HermiteCoulombIntegral", "An implementation of the (auxiliary) Hermite Coulomb integral R^n_{tuv}."};

    py_HermiteCoulombIntegral
        .def(py::init<const double, const Vector<double, 3>, const Vector<double, 3>>(),
             py::arg("p"),
             py::arg("P"),
             py::arg("C"))

        .def("__call__",
             &HermiteCoulombIntegral::operator(),
             py::arg("n"),
             py::arg("t"),
             py::arg("u"),
             py::arg("v"),
             "Return the value for the (auxiliary) Hermite Coulomb integral R^n_{tuv}(p, P, C).");
}


}  // namespace gqcpy
