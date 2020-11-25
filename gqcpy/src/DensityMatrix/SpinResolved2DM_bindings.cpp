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

#include "DensityMatrix/SpinResolved2DM.hpp"
#include "gqcpy/include/interfaces.hpp"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>


namespace py = pybind11;


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


/**
 *  Register `SpinResolved2DM_d` to the gqcpy module and expose a part of its C++ interface to Python.
 * 
 *  @param module           The Pybind11 module in which `SpinResolved2DM_d` should be registered.
 */
void bindSpinResolved2DM(py::module& module) {

    // Define the Python class for `SpinResolved2DM`.
    py::class_<SpinResolved2DM<double>> py_SpinResolved2DM_d {module, "SpinResolved2DM_d", "A type that encapsulates the spin parts of the spin-resolved two-electron density matrix."};


    py_SpinResolved2DM_d

        /*
         *  MARK: General information
         */

        .def(
            "numberOfOrbitals",
            &SpinResolved2DM<double>::numberOfOrbitals,
            "Return the number of orbitals (spinors or spin-orbitals, depending on the context) that are related to the sigma-tau part of the spin-resolved 2-DM.")

        .def(
            "orbitalDensity",
            &SpinResolved2DM<double>::orbitalDensity,
            "The orbital (total, spin-summed) two-electron density matrix.");


    // Expose the `DoublySpinResolvedBase` API to the Python class.
    bindDoublySpinResolvedBaseInterface(py_SpinResolved2DM_d);
}


}  // namespace gqcpy
