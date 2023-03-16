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

#include "Domain/MullikenDomain/GMullikenDomain.hpp"
#include "Utilities/complex.hpp"
#include "gqcpy/include/interfaces.hpp"

#include <pybind11/pybind11.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


/**
 *  Register `GMullikenDomain` to the gqcpy module and expose a part of its C++ interface to Python.
 *
 *  @param module           The Pybind11 module in which `GMullikenDomain` should be registered.
 */
void bindGMullikenDomain(py::module& module) {

    // Define the Python class for `GMullikenDomain_d`.
    py::class_<GMullikenDomain<double>> py_GMullikenDomain_d {module, "GMullikenDomain_d", "A real generalized Mulliken-based domain of an AO basis."};

    // Expose the "Mulliken Domain" interface to the Python class.
    bindMullikenDomainIndicesInterface(py_GMullikenDomain_d);
    bindMullikenDomainMatricesInterface(py_GMullikenDomain_d);

    // Define the Python class for `GMullikenDomain_cd`.
    py::class_<GMullikenDomain<complex>> py_GMullikenDomain_cd {module, "GMullikenDomain_cd", "A complex generalized Mulliken-based domain of an AO basis."};

    // Expose the "Mulliken Domain" interface to the Python class.
    bindMullikenDomainIndicesInterface(py_GMullikenDomain_cd);
    bindMullikenDomainMatricesInterface(py_GMullikenDomain_cd);
}


}  // namespace gqcpy
