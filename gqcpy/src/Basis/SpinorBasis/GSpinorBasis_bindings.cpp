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

#include "Basis/ScalarBasis/GTOShell.hpp"
#include "Basis/SpinorBasis/GSpinorBasis.hpp"
#include "Molecule/Molecule.hpp"
#include "Operator/FirstQuantized/Operator.hpp"
#include "gqcpy/include/interfaces.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>


namespace py = pybind11;


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


/**
 *  Register `GSpinorBasis_d` to the gqcpy module and expose a part of its C++ interface to Python.
 * 
 *  @param module           The Pybind11 module in which `GSpinorBasis_d` should be registered.
 */
void bindGSpinorBasis(py::module& module) {

    // Define the Python class for `GSpinorBasis`.
    py::class_<GSpinorBasis<double, GTOShell>> py_GSpinorBasis_d {module, "GSpinorBasis_d", "A class that represents a real, (generalized) spinor basis with underlying GTO shells."};


    /*
     *  MARK: Constructors
     */

    py_GSpinorBasis_d
        .def(py::init<const Molecule&, const std::string&>(),
             py::arg("molecule"),
             py::arg("basisset_name"));


    /*
     *  MARK: Named constructors
     */

    py_GSpinorBasis_d
        .def_static("FromRestricted",
                    [](const RSpinOrbitalBasis<double, GTOShell>& r_spinor_basis) {
                        return GSpinorBasis<double, GTOShell>::FromRestricted(r_spinor_basis);
                    });


    /*
     *  MARK: Coefficients
     */

    py_GSpinorBasis_d
        .def(
            "numberOfCoefficients",
            &GSpinorBasis<double, GTOShell>::numberOfCoefficients,
            py::arg("sigma"),
            "Return the number of coefficients that are used for the expansion of the requested spin-component of a spinor");


    // Expose the `SimpleSpinorBasis` API to the Python class.
    bindSpinorBasisInterface(py_GSpinorBasis_d);


    /*
     *  MARK: General info
     */

    py_GSpinorBasis_d
        .def(
            "numberOfSpinors",
            &GSpinorBasis<double, GTOShell>::numberOfSpinors,
            "The number of spinors that are described by this generalized spinor basis.");


    // Expose some quantization API to the Python class;
    bindSpinorBasisQuantizationInterface(py_GSpinorBasis_d);


    // Expose some Mulliken API to the Python class;
    bindSpinorBasisMullikenInterface(py_GSpinorBasis_d);
}


}  // namespace gqcpy
