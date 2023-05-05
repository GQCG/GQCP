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

#include "Partition/SpinResolvedElectronPartition.hpp"
#include "gqcpy/include/interfaces.hpp"

#include <pybind11/functional.h>
#include <pybind11/pybind11.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


void bindSpinResolvedElectronPartition(py::module& module) {
    py::class_<SpinResolvedElectronPartition> py_SpinResolvedElectronPartition {module, "SpinResolvedElectronPartition", "A spin-resolved partition of alpha and beta electron numbers over e.g. domains."};


    py_SpinResolvedElectronPartition

        /*
         *  MARK: General info
         */
        .def(
            "__eq__",
            [](const SpinResolvedElectronPartition& electron_partition, const SpinUnresolvedElectronPartition& other) {
                return electron_partition == other;
            },
            py::arg("other"),
            "Return whether this `SpinResolvedElectronPartition` is equal to a `SpinUnresolvedElectronPartition`.")

        .def(
            "__eq__",
            [](const SpinResolvedElectronPartition& electron_partition, const SpinResolvedElectronPartition& other) {
                return electron_partition == other;
            },
            py::arg("other"),
            "Return whether this `SpinResolvedElectronPartition` is equal to a `SpinResolvedElectronPartition`.")

        .def(
            "asString",
            &SpinResolvedElectronPartition::asString,
            "Return the spin-resolved electron partition string representation.")

        .def(
            "numberOfElectrons",
            [](const SpinResolvedElectronPartition& electron_partition, size_t i) {
                return electron_partition.numberOfElectrons(i);
            },
            py::arg("i"),
            "Return the number of alpha and beta electrons the partition contains at index `i`.")

        .def(
            "numberOfElectrons",
            [](const SpinResolvedElectronPartition& electron_partition) {
                return electron_partition.numberOfElectrons();
            },
            "Return the total number of electrons the partition contains.");


    // Expose the `SpinResolvedBase` interface to the Python class.
    bindSpinResolvedBaseInterface(py_SpinResolvedElectronPartition);
}


}  // namespace gqcpy
