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

#include "DensityMatrix/SpinResolved1DM.hpp"
#include "gqcpy/include/interfaces.hpp"

#include <pybind11/pybind11.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


/**
 *  Register `SpinResolved1DM_d` to the gqcpy module and expose a part of its C++ interface to Python.
 * 
 *  @param module           The Pybind11 module in which `SpinResolved1DM_d` should be registered.
 */
void bindSpinResolved1DM(py::module& module) {

    // Define the Python class for `SpinResolved1DM`.
    py::class_<SpinResolved1DM<double>> py_SpinResolved1DM_d {module, "SpinResolved1DM_d", "A type that encapsulates alpha and beta (spin-resolved) density matrices."};


    py_SpinResolved1DM_d

        /*
         *  MARK: Named constructors
         */

        .def_static(
            "FromOrbital1DM",
            &SpinResolved1DM<double>::FromOrbital1DM,
            "Create a spin-resolved 1-DM from an `Orbital1DM`, attributing half of the orbital 1-DM to each of the spin components.")


        /*
         *  MARK: General information
         */

        .def(
            "numberOfOrbitals",
            &SpinResolved1DM<double>::numberOfOrbitals,
            py::arg("sigma"),
            "Return the number of orbitals (spinors or spin-orbitals, depending on the context) that correspond to the given spin.")


        /*
         *  MARK: Spin-related operations
         */

        .def(
            "spinDensity",
            &SpinResolved1DM<double>::spinDensity,
            "Return the spin-density matrix, i.e. the difference between the alpha and beta 1-DM.")

        .def(
            "orbitalDensity",
            &SpinResolved1DM<double>::orbitalDensity,
            "Return the orbital density matrix, i.e. the sum of the alpha and beta 1-DM.");


    // Expose the `SpinResolvedBase` API to the Python class.
    bindSpinResolvedBaseInterface(py_SpinResolved1DM_d);


    // Expose the `BasisTransformable` API to the Python class.
    bindBasisTransformableInterface(py_SpinResolved1DM_d);

    // Expose the `VectorSpaceArithmetic` API to the Python class.
    bindVectorSpaceArithmeticInterface(py_SpinResolved1DM_d);
}


}  // namespace gqcpy
