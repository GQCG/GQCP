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

#include "ONVBasis/SpinUnresolvedONVBasis.hpp"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


void bindSpinUnresolvedONVBasis(py::module& module) {
    py::class_<SpinUnresolvedONVBasis>(module, "SpinUnresolvedONVBasis", "A full spin-unresolved ONV basis.")

        // CONSTRUCTORS

        .def(py::init<const size_t, const size_t>(),
             py::arg("M"),
             py::arg("N"))


        // PUBLIC METHODS

        .def(
            "arcWeight",
            [](const SpinUnresolvedONVBasis& onv_basis, size_t p, const size_t n) {
                return onv_basis.arcWeight(p, n);
            },
            py::arg("p"),
            py::arg("n"),
            "Return the arc weight of the arc starting at a vertex (p, n), with p the orbital index and n the electron index.")

        .def(
            "calculateDimension",
            [](const SpinUnresolvedONVBasis& onv_basis, size_t M, const size_t N) {
                return onv_basis.calculateDimension(M, N);
            },
            py::arg("M"),
            py::arg("N"),
            "Calculate the dimension of a spin-unresolved ONV with M spinors and N electrons.")

        .def(
            "vertexWeight",
            [](const SpinUnresolvedONVBasis onv_basis, const size_t p, const size_t n) {
                return onv_basis.vertexWeight(p, n);
            },
            py::arg("p"),
            py::arg("n"),
            "Return the vertex weight of vertex (p, n), with p the orbital index and n the electron index.");
}

}  // namespace gqcpy
