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

#include "Molecule/Nucleus.hpp"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


void bindNucleus(py::module& module) {
    py::class_<Nucleus>(module, "Nucleus", "A class that represents a nucleus: it has a charge and a position in space")

        /*
         *  MARK: Constructors
         */

        .def(py::init<const size_t, const double, const double, const double>(),
             py::arg("Z"),
             py::arg("x"),
             py::arg("y"),
             py::arg("z"))


        /*
         *  MARK: General information
         */

        .def("element",
             &Nucleus::element,
             "Return the string representation of the element that corresponds to this nucleus.")

        .def("position",
             &Nucleus::position,
             "Return the vector describing the position of the nucleus.")

        .def("charge",
             &Nucleus::charge,
             "Return the nuclear charge");
}


}  // namespace gqcpy
