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

#include "Domain/MullikenDomain/UMullikenDomainComponent.hpp"
#include "Utilities/complex.hpp"
#include "gqcpy/include/interfaces.hpp"

#include <pybind11/pybind11.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


/**
 *  Register `UMullikenDomainComponent` to the gqcpy module and expose a part of its C++ interface to Python.
 *
 *  @param module           The Pybind11 module in which `UMullikenDomainComponent` should be registered.
 */
void bindUMullikenDomainComponent(py::module& module) {

    // Define the Python class for `UMullikenDomainComponent_d`.
    py::class_<UMullikenDomainComponent<double>> py_UMullikenDomainComponent_d {module, "UMullikenDomainComponent_d", "One of the components of a real unrestricted Mulliken-based domain of an AO basis."};

    // Expose the "Mulliken Domain" interface to the Python class.
    bindDiscreteDomainInterface(py_UMullikenDomainComponent_d);
    bindMullikenDomainMatricesInterface(py_UMullikenDomainComponent_d);

    // Define the Python class for `UMullikenDomainComponent_cd`.
    py::class_<UMullikenDomainComponent<complex>> py_UMullikenDomainComponent_cd {module, "UMullikenDomainComponent_cd", "One of the components of a complex unrestricted Mulliken-based domain of an AO basis."};

    // Expose the "Mulliken Domain" interface to the Python class.
    bindDiscreteDomainInterface(py_UMullikenDomainComponent_cd);
    bindMullikenDomainMatricesInterface(py_UMullikenDomainComponent_cd);
}


}  // namespace gqcpy
