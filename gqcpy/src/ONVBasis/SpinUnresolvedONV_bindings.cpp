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

#include "ONVBasis/SpinUnresolvedONV.hpp"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


void bindSpinUnresolvedONV(py::module& module) {
    py::class_<SpinUnresolvedONV>(module, "SpinUnresolvedONV", "A spin-unresolved occupation number vector.")

        // CONSTRUCTORS

        .def_static(
            "FromString",
            [](const std::string& string_representation) {
                return SpinUnresolvedONV::FromString(string_representation);
            },
            py::arg("string_representation"),
            "Create a spin-unresolved ONV from a string representation.")

        .def_static(
            "FromOccupiedIndices",
            [](const std::vector<size_t>& occupied_indices, const size_t M) {
                return SpinUnresolvedONV::FromOccupiedIndices(occupied_indices, M);
            },
            py::arg("occupied_indices"),
            py::arg("M"),
            "Create a spin-unresolved ONV from an array of occupied indices and a total number of spinors.")

        .def_static(
            "GHF",
            [](const size_t M, const size_t N, const Eigen::MatrixXd& orbital_energies) {
                return SpinUnresolvedONV::GHF(M, N, VectorX<double>(orbital_energies));
            },
            py::arg("M"),
            py::arg("N"),
            py::arg("orbital_energies"),
            "Create a spin-unresolved ONV that represents the GHF single Slater determinant.")


        // PUBLIC METHODS

        .def(
            "__repr__",
            &SpinUnresolvedONV::asString)

        .def(
            "calculateProjection",
            [](const SpinUnresolvedONV& onv_of, const SpinUnresolvedONV& onv_on, const Eigen::MatrixXd& C_of, const Eigen::MatrixXd& C_on, const Eigen::MatrixXd& S) {
                return onv_of.calculateProjection(onv_on, GTransformation<double>(C_of), GTransformation<double>(C_on), SquareMatrix<double>(S));
            },
            py::arg("onv_on"),
            py::arg("C_of"),
            py::arg("C_on"),
            py::arg("S"),
            "Calculate the overlap <on|of>: the projection of between this spin-unresolved ONV ('of') and another spin-unresolved ONV ('on'), expressed in different general orthonormal spinor bases.")

        .def(
            "orbitalSpace",
            &SpinUnresolvedONV::orbitalSpace,
            "Return the implicit orbital space that is related to this spin-unresolved ONV by taking this as a reference determinant.");
}


}  // namespace gqcpy
