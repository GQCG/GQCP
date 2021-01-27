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


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


void bindHoppingMatrix(py::module& module) {
    py::class_<HoppingMatrix<double>>(module, "HoppingMatrix", "The Hubbard hopping matrix.")

        /*
         *  MARK: Constructors
         */

        .def(py::init<const AdjacencyMatrix&, const double, const double>(),
             py::arg("A"),
             py::arg("t"),
             py::arg("U"),
             "Return the Hubbard hopping matrix from an adjacency matrix and Hubbard model parameters U and t.")

        .def_static(
            "FromCSLine",
            [](const std::string& cs_line) {
                return HoppingMatrix<double>::FromCSLine(cs_line);
            },
            py::arg("csline"),
            "Return the hopping matrix that corresponds to the given comma-separated line.")

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
