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

#include "ONVBasis/SeniorityZeroONVBasis.hpp"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


void bindSeniorityZeroONVBasis(py::module& module) {
    py::class_<SeniorityZeroONVBasis>(module, "SeniorityZeroONVBasis", "A full seniority-zero ONV basis.")

        // CONSTRUCTORS

        .def(py::init<const size_t, const size_t>(),
             py::arg("K"),
             py::arg("N_P"))

        // PUBLIC METHODS

        .def("dimension",
             &SeniorityZeroONVBasis::dimension);
}


}  // namespace gqcpy
