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

#include "Domain/HubbardDomain.hpp"
#include "ONVBasis/SpinResolvedONV.hpp"
#include "ONVBasis/SpinUnresolvedONV.hpp"
#include "Utilities/complex.hpp"
#include "gqcpy/include/interfaces.hpp"

#include <pybind11/pybind11.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


/**
 *  Register `HubbardDomain` to the gqcpy module and expose a part of its C++ interface to Python.
 *
 *  @param module           The Pybind11 module in which `HubbardDomain` should be registered.
 */
void bindHubbardDomain(py::module& module) {

    // Define the Python class for `HubbardDomain`.
    py::class_<HubbardDomain> py_HubbardDomain {module, "HubbardDomain", "A real Hubbard domain as a partial set of the sites of the system."};

    // Expose python bindings unique to the Hubbard domain.
    py_HubbardDomain

        .def(
            "overlapWithONV",
            [](const HubbardDomain& domain, const GQCP::SpinUnresolvedONV& onv) {
                return domain.overlapWithONV(onv);
            },
            py::arg("unresolved_ONV"),
            "Return the overlap of this domain with a given spin-unresolved ONV.")

        .def(
            "overlapWithONV",
            [](const HubbardDomain& domain, const GQCP::SpinResolvedONV& onv) {
                return domain.overlapWithONV(onv);
            },
            py::arg("resolved_ONV"),
            "Return the overlap of this domain with a given spin-resolved ONV.");


    // Expose the "Mulliken Domain" interface to the Python class.
    bindDiscreteDomainInterface(py_HubbardDomain);
}


}  // namespace gqcpy
