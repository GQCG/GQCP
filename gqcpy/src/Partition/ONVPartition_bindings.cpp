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

#include "Partition/ONVPartition.hpp"
#include "gqcpy/include/interfaces.hpp"

#include <pybind11/pybind11.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


/**
 *  Register `ONVPartition` to the gqcpy module and expose a part of its C++ interface to Python.
 *
 *  @tparam ONV               The type of ONV.
 *
 *  @param py_class             The Pybind11 `class_` that should obtain APIs related to `ONVPartition`.
 */
template <typename ONV>
void bindONVPartitionInterface(py::module& module, const std::string& suffix, const std::string& description) {

    // Define the Python class for `DiscreteDomainPartition`.
    py::class_<ONVPartition<ONV>> py_ONVPartition(module, ("ONVPartition_" + suffix).c_str(), description.c_str());

    // Expose python bindings unique to the ONV partition.
    py_ONVPartition

        /*
         *  MARK: Constructors
         */

        .def(py::init<const DiscreteDomainPartition&, const ONV&>(),
             py::arg("domain_partition"),
             py::arg("onv"))

        /*
         *  MARK: General info
         */

        .def(
            "asString",
            &ONVPartition<ONV>::asString)

        .def(
            "phaseFactor",
            &ONVPartition<ONV>::phaseFactor);


    // Expose the `SimplePartition` interfaces to the Python class.
    bindSimplePartitionInterface(py_ONVPartition);
}


/**
 *  Register multiple templated `ONVPartition` to the gqcpy module and expose parts of their C++ interfaces to Python.
 *
 *  @param module           The Pybind11 module in which the classes should be registered.
 */
void bindONVPartition(py::module& module) {

    bindONVPartitionInterface<SpinResolvedONV>(module, "SpinResolved", "A spin-resolved ONV partitioned according to a set of discrete domains, i.e. a `DiscreteDomainPartition`.");
    bindONVPartitionInterface<SpinUnresolvedONV>(module, "SpinUnresolved", "A spin-unresolved ONV partitioned according to a set of discrete domains, i.e. a `DiscreteDomainPartition`.");
}


}  // namespace gqcpy
