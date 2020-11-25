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
#include "Basis/SpinorBasis/RSpinOrbitalBasis.hpp"
#include "Operator/FirstQuantized/Operator.hpp"
#include "gqcpy/include/interfaces.hpp"

#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


/**
 *  Register `RSpinOrbitalBasis_d` to the gqcpy module and expose a part of its C++ interface to Python.
 * 
 *  @param module           The Pybind11 module in which `RSpinOrbitalBasis_d` should be registered.
 */
void bindRSpinOrbitalBasis(py::module& module) {

    // Define the Python class for `RSpinOrbitalBasis`.
    py::class_<RSpinOrbitalBasis<double, GTOShell>> py_RSpinOrbitalBasis_d {module, "RSpinOrbitalBasis_d", "A class that represents a real, restricted spinor basis with underlying GTO shells."};

    /*
     *  MARK: Constructors
     */

    py_RSpinOrbitalBasis_d
        .def(py::init<const Molecule&, const std::string&>(),
             py::arg("molecule"),
             py::arg("basisset_name"));


    // Expose the `SimpleSpinorBasis` API to the Python class.
    bindSpinorBasisInterface(py_RSpinOrbitalBasis_d);

    /*
     *  MARK: General information
     */

    py_RSpinOrbitalBasis_d
        .def(
            "numberOfSpatialOrbitals",
            [](RSpinOrbitalBasis<double, GTOShell>& spin_orbital_basis) {
                return spin_orbital_basis.numberOfSpatialOrbitals();
            },
            "Return the number of different spatial orbitals that are used in this spinor basis.");


    /*
     *  MARK: Quantization of first-quantized operators
     */

    py_RSpinOrbitalBasis_d
        .def(
            "quantizeDipoleOperator",
            [](const RSpinOrbitalBasis<double, GTOShell>& spin_orbital_basis, const Vector<double, 3>& origin) {
                return spin_orbital_basis.quantize(Operator::ElectronicDipole(origin));
            },
            py::arg("origin") = Vector<double, 3>::Zero(), "Return the electronic dipole operator expressed in this spinor basis.");


    // Expose some quantization API to the Python class;
    bindSpinorBasisQuantizationInterface(py_RSpinOrbitalBasis_d);

    // Expose some Mulliken API to the Python class;
    bindSpinorBasisMullikenInterface(py_RSpinOrbitalBasis_d);
}


}  // namespace gqcpy
