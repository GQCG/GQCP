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

#include "Domain/DiscreteDomain.hpp"
#include "ONVBasis/SpinResolvedONV.hpp"
#include "ONVBasis/SpinUnresolvedONV.hpp"
#include "gqcpy/include/interfaces.hpp"

#include <pybind11/pybind11.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


/**
 *  Register `DiscreteDomain` to the gqcpy module and expose a part of its C++ interface to Python.
 *
 *  @param module           The Pybind11 module in which `DiscreteDomain` should be registered.
 */
void bindDiscreteDomain(py::module& module) {

    // Define the Python class for `DiscreteDomain`.
    py::class_<DiscreteDomain> py_DiscreteDomain {module, "DiscreteDomain", "A domain represented by a bitstring that specifies whether the element at index `i' is present in the domain or not."};

    // Expose python bindings unique to the discrete domain.
    py_DiscreteDomain

        /*
         * MARK: Overlap
         */

        .def(
            "overlapWith",
            [](const DiscreteDomain& domain, const DiscreteDomain& other_domain) {
                return domain.overlapWith(other_domain);
            },
            py::arg("other_domain"),
            "Return the number of domain elements that are equal between this domain and the other discrete domain.")

        .def(
            "overlapWithONV",
            [](const DiscreteDomain& domain, const GQCP::SpinUnresolvedONV& onv) {
                return domain.overlapWithONV(onv);
            },
            py::arg("unresolved_ONV"),
            "Return the overlap of this domain with a given spin-unresolved ONV.")

        .def(
            "overlapWithONV",
            [](const DiscreteDomain& domain, const GQCP::SpinResolvedONV& onv) {
                return domain.overlapWithONV(onv);
            },
            py::arg("resolved_ONV"),
            "Return the overlap of this domain with a given spin-resolved ONV.");


    // Expose the "DiscreteDomain" interfaces to the Python class.
    bindDiscreteDomainConstructorInterface(py_DiscreteDomain);
    bindDiscreteDomainInterface(py_DiscreteDomain);
}


}  // namespace gqcpy
