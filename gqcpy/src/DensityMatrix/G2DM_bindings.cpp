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

#include "DensityMatrix/G2DM.hpp"
#include "gqcpy/include/interfaces.hpp"

#include <pybind11/pybind11.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


/**
 *  Register `G2DM_d` to the gqcpy module and expose a part of its C++ interface to Python.
 * 
 *  @param module           The Pybind11 module in which `G2DM_d` should be registered.
 */
void bindG2DM(py::module& module) {

    // Define the Python class for `G2DM`.
    py::class_<G2DM<double>> py_G2DM_d {module, "G2DM_d", "A type used to represent a two-electron general(ized) density matrix, i.e. the full spinor two-component two-electron density matrix."};


    // Expose the `Simple2DM` API to the Python class.
    bindSimple2DMInterface(py_G2DM_d);
}


}  // namespace gqcpy
