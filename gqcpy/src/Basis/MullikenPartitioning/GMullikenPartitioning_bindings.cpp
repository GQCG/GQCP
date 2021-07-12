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

#include "Basis/MullikenPartitioning/GMullikenPartitioning.hpp"
#include "gqcpy/include/interfaces.hpp"
#include "Utilities/complex.hpp"

#include <pybind11/pybind11.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


/**
 *  Register `GMullikenPartitioning` to the gqcpy module and expose a part of its C++ interface to Python.
 * 
 *  @param module           The Pybind11 module in which `GMullikenPartitioning` should be registered.
 */
void bindGMullikenPartitioning(py::module& module) {

    // Define the Python class for `GMullikenPartitioning_d`.
    py::class_<GMullikenPartitioning<double>> py_GMullikenPartitioning_d {module, "GMullikenPartitioning_d", "A real generalized Mulliken-based partitioning of an AO basis."};

    // Expose the "Mulliken partitioning" interface to the Python class.
    bindMullikenPartitioningIndicesInterface(py_GMullikenPartitioning_d);
    bindMullikenPartitioningMatricesInterface(py_GMullikenPartitioning_d);

    // Define the Python class for `GMullikenPartitioning_cd`.
    py::class_<GMullikenPartitioning<complex>> py_GMullikenPartitioning_cd {module, "GMullikenPartitioning_cd", "A complex generalized Mulliken-based partitioning of an AO basis."};

    // Expose the "Mulliken partitioning" interface to the Python class.
    bindMullikenPartitioningIndicesInterface(py_GMullikenPartitioning_cd);
    bindMullikenPartitioningMatricesInterface(py_GMullikenPartitioning_cd);
}


}  // namespace gqcpy
