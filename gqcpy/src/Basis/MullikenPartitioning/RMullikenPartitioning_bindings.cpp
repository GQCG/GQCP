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

#include "Basis/MullikenPartitioning/RMullikenPartitioning.hpp"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>


namespace py = pybind11;


namespace gqcpy {

using namespace GQCP;


void bindRMullikenPartitioning(py::module& module) {
    py::class_<RMullikenPartitioning<double>>(module, "RMullikenPartitioning", "A restricted Mulliken-based partitioning of an AO basis.")

        /*
         *  MARK: General information
         */

        .def("indices",
             &RMullikenPartitioning<double>::indices,
             "Return the set of indices that correspond to the AOs that are included in the Mulliken-partitioning of an AO basis.")

        /*
         *  MARK: Partitioning and projecting
         */

        .def("partitionMatrix",
             &RMullikenPartitioning<double>::partitionMatrix,
             "Return the partition matrix 'P_A' related to this Mulliken partitioning.")

        .def("projectionMatrix",
             &RMullikenPartitioning<double>::projectionMatrix,
             "Return the Mulliken projection, defined as C^{-1} P_A C, where C is the transformation matrix and P_A is the partition matrix.");
}


}  // namespace gqcpy
