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

#include "Partition/SpinUnresolvedElectronPartition.hpp"
#include "gqcpy/include/interfaces.hpp"

#include <pybind11/pybind11.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


/**
 *  Register `SpinUnresolvedElectronPartition` to the gqcpy module and expose a part of its C++ interface to Python.
 *
 *  @param module           The Pybind11 module in which `SpinUnresolvedElectronPartition` should be registered.
 */
void bindSpinUnresolvedElectronPartition(py::module& module) {

    // Define the Python class for `SpinUnresolvedElectronPartition`.
    py::class_<SpinUnresolvedElectronPartition> py_SpinUnresolvedElectronPartition {module, "SpinUnresolvedElectronPartition", "A partition of electron numbers over e.g. domains."};

    // Expose python bindings unique to the spin-unresolved electron partition.
    py_SpinUnresolvedElectronPartition

        /*
         *  MARK: General info
         */

        .def(
            "asString",
            &SpinUnresolvedElectronPartition::asString)

        .def(
            "numberOfElectrons",
            [](const SpinUnresolvedElectronPartition& electron_partition, size_t i) {
                return electron_partition.numberOfElectrons(i);
            },
            py::arg("i"),
            "Return the number of electrons the partition contains at index `i`.")

        .def(
            "numberOfElectrons",
            [](const SpinUnresolvedElectronPartition& electron_partition) {
                return electron_partition.numberOfElectrons();
            },
            "Return the total number of electrons the partition contains.");


    // Expose the `SimplePartition` interfaces to the Python class.
    bindSimplePartitionInterface(py_SpinUnresolvedElectronPartition);
}


}  // namespace gqcpy
