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

#include "Domain/MullikenDomain/RMullikenDomain.hpp"
#include "Utilities/complex.hpp"
#include "gqcpy/include/interfaces.hpp"

#include <pybind11/pybind11.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


/**
 *  Register `RMullikenDomain` to the gqcpy module and expose a part of its C++ interface to Python.
 *
 *  @param module           The Pybind11 module in which `RMullikenDomain` should be registered.
 */
void bindRMullikenDomain(py::module& module) {

    // Define the Python class for `RMullikenDomain_d`.
    py::class_<RMullikenDomain<double>> py_RMullikenDomain_d {module, "RMullikenDomain_d", "A real restricted Mulliken-based domain of an AO basis."};

    // Expose the "Mulliken Domain" interface to the Python class.
    bindDiscreteDomainInterface(py_RMullikenDomain_d);
    bindMullikenDomainMatricesInterface(py_RMullikenDomain_d);

    // Define the Python class for `RMullikenDomain_cd`.
    py::class_<RMullikenDomain<complex>> py_RMullikenDomain_cd {module, "RMullikenDomain_cd", "A complex restricted Mulliken-based domain of an AO basis."};

    // Expose the "Mulliken Domain" interface to the Python class.
    bindDiscreteDomainInterface(py_RMullikenDomain_cd);
    bindMullikenDomainMatricesInterface(py_RMullikenDomain_cd);
}


}  // namespace gqcpy
