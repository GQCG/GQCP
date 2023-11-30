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

#include "Operator/SecondQuantized/ModelHamiltonian/HoppingMatrix.hpp"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


void bindHoppingMatrix(py::module& module) {
    py::class_<HoppingMatrix<double>>(module, "HoppingMatrix", "The Hubbard hopping matrix.")

        /*
         *  MARK: Constructors
         */

        .def_static(
            "Homogeneous",
            [](const AdjacencyMatrix& A, const double t) {
                return HoppingMatrix<double>::Homogeneous(A, t);
            },
            py::arg("A"),
            py::arg("t"),

            "Return the Hubbard hopping matrix from an adjacency matrix and Hubbard model parameter t.")

        .def_static(
            "Dense",
            [](std::vector<double>& triagonal_data) {
                return HoppingMatrix<double>::Dense(triagonal_data);
            },
            py::arg("triagonal_data"),
            "Return the hopping matrix that corresponds to the given upper triangle values.")

        .def_static(
            "FromLinkVector",
            [](const AdjacencyMatrix& A, std::vector<double>& link_vector) {
                return HoppingMatrix<double>::FromLinkVector(A, link_vector);
            },
            py::arg("A"),
            py::arg("link_vector"),
            "Return the Hubbard hopping matrix from an adjacency matrix and link vector containing values for different parameters t.")

        /*
         *  MARK: Access
         */
        .def(
            "matrix",
            [](const HoppingMatrix<double>& H) {
                return H.matrix().Eigen();
            },
            "Return a read-only reference to the matrix representation of this hopping matrix.")

        /*
         *  MARK: General information
         */
        .def(
            "numberOfLatticeSites",
            [](const HoppingMatrix<double>& H) {
                return H.numberOfLatticeSites();
            },
            "Return the number of lattice sites corresponding used in this hopping matrix.");
}


}  // namespace gqcpy
