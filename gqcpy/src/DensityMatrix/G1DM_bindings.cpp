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

#include "DensityMatrix/G1DM.hpp"
#include "gqcpy/include/interfaces.hpp"
#include "Utilities/aliases.hpp"

#include <pybind11/pybind11.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


/**
 *  Register `G1DM_d` and `G1DM_cd` to the gqcpy module and expose a part of their C++ interface to Python.
 * 
 *  @param module           The Pybind11 module in which `G1DM_d` and `G2DM_cd` should be registered.
 */
void bindG1DM(py::module& module) {

    // Define the Python class for `G1DM_d`.
    py::class_<G1DM<double>> py_G1DM_d {module, "G1DM_d", "A type used to represent a real-valued one-electron general(ized) density matrix, i.e. the full spinor two-component one-electron density matrix."};

    // Expose the `Simple1DM` and `BasisTransformable` APIs to `G1DM_d;
    bindSimple1DMInterface(py_G1DM_d);
    bindBasisTransformableInterface(py_G1DM_d);
    

    // Define the Python class for `G1DM_cd`.
    py::class_<G1DM<complex>> py_G1DM_cd {module, "G1DM_cd", "A type used to represent a complex-valued one-electron general(ized) density matrix, i.e. the full spinor two-component one-electron density matrix."};

    // Expose the `Simple1DM` and `BasisTransformable` APIs to `G1DM_cd;
    bindSimple1DMInterface(py_G1DM_cd);
    bindBasisTransformableInterface(py_G1DM_cd);
}


}  // namespace gqcpy
