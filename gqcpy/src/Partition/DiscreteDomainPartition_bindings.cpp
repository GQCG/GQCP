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

#include "Partition/DiscreteDomainPartition.hpp"
#include "gqcpy/include/interfaces.hpp"

#include <pybind11/pybind11.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


/**
 *  Register `DiscreteDomainPartition` to the gqcpy module and expose a part of its C++ interface to Python.
 *
 *  @param module           The Pybind11 module in which `DiscreteDomainPartition` should be registered.
 */
void bindDiscreteDomainPartition(py::module& module) {

    // Define the Python class for `DiscreteDomainPartition`.
    py::class_<DiscreteDomainPartition> py_DiscreteDomainPartition {module, "DiscreteDomainPartition", "A partition (i.e., collection) of discrete domains."};

    // Expose python bindings unique to the discrete domain partition.
    py_DiscreteDomainPartition

        /*
         *  MARK: General info
         */
        .def(
            "asString",
            &DiscreteDomainPartition::asString)

        .def(
            "asVector",
            &DiscreteDomainPartition::asVector)

        /*
         * MARK: Overlap
         */

        .def(
            "overlapWithONV",
            [](const DiscreteDomainPartition& domain_partition, const GQCP::SpinUnresolvedONV& onv) {
                return domain_partition.overlapWithONV(onv);
            },
            py::arg("unresolved_ONV"),
            "Return the numbers of overlapping set bits after a bit-by-bit comparison between the discrete domains and the spin-unresolved ONV.")

        .def(
            "overlapWithONV",
            [](const DiscreteDomainPartition& domain_partition, const GQCP::SpinResolvedONV& onv) {
                return domain_partition.overlapWithONV(onv);
            },
            py::arg("resolved_ONV"),
            "Return the numbers of overlapping set bits after a bit-by-bit comparison between the discrete domains and the spin-resolved ONV.");


    // Expose the `DomainPartition` interfaces to the Python class.
    bindDomainPartitionInterface(py_DiscreteDomainPartition);
}


}  // namespace gqcpy
