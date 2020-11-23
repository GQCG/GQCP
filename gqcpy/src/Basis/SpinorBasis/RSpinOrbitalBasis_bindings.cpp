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

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>


namespace py = pybind11;


namespace gqcpy {

using namespace GQCP;


void bindRSpinOrbitalBasis(py::module& module) {
    py::class_<RSpinOrbitalBasis<double, GTOShell>>(module, "RSpinOrbitalBasis", "A class that represents a real, restricted spinor basis with underlying GTO shells.")

        /*
         *  MARK: Constructors
         */

        .def(py::init<const Molecule&, const std::string&>(),
             py::arg("molecule"),
             py::arg("basisset_name"))


        /*
         *  MARK: Inherited methods
         */

        .def(
            "expansion",
            &RSpinOrbitalBasis<double, GTOShell>::expansion,
            "Return a read-only reference to the transformation that relates the current set of spinors with the atomic spinors.")

        .def(
            "overlap",
            &RSpinOrbitalBasis<double, GTOShell>::overlap,
            "Return the overlap (one-electron) operator of this spin-orbital basis.")

        .def(
            "isOrthonormal",
            &RSpinOrbitalBasis<double, GTOShell>::isOrthonormal,
            py::arg("precision") = 1.0e-08,
            "If this spinor basis is orthonormal.")

        .def(
            "lowdinOrthonormalization",
            &RSpinOrbitalBasis<double, GTOShell>::lowdinOrthonormalization,
            "Return the transformation to the Löwdin basis: T = S_current^{-1/2}.")

        .def(
            "lowdinOrthonormalize",
            &RSpinOrbitalBasis<double, GTOShell>::lowdinOrthonormalize,
            " Transform the spinor basis to the 'Löwdin basis', which is the orthonormal basis that we transform to with T = S^{-1/2}, where S is the current overlap matrix.")


        /*
         *  MARK: General information
         */

        .def(
            "numberOfSpatialOrbitals",
            [](RSpinOrbitalBasis<double, GTOShell>& spin_orbital_basis) {
                return spin_orbital_basis.numberOfSpatialOrbitals();
            },
            "Return the number of different spatial orbitals that are used in this spinor basis.")


        /*
         *  MARK: Quantization of first-quantized operators
         */

        .def(
            "quantizeDipoleOperator",
            [](const RSpinOrbitalBasis<double, GTOShell>& spin_orbital_basis, const Vector<double, 3>& origin) {
                return spin_orbital_basis.quantize(Operator::ElectronicDipole(origin));
            },
            py::arg("origin") = Vector<double, 3>::Zero(), "Return the electronic dipole operator expressed in this spinor basis.")

        .def(
            "quantizeKineticOperator",
            [](const RSpinOrbitalBasis<double, GTOShell>& spin_orbital_basis) {
                return spin_orbital_basis.quantize(Operator::Kinetic());
            },
            "Return the kinetic energy operator expressed in this spinor basis.")

        .def(
            "quantizeNuclearAttractionOperator",
            [](const RSpinOrbitalBasis<double, GTOShell>& spin_orbital_basis, const Molecule& molecule) {
                return spin_orbital_basis.quantize(Operator::NuclearAttraction(molecule));
            },
            "Return the nuclear attraction operator expressed in this spinor basis.")

        .def(
            "quantizeOverlapOperator",
            [](const RSpinOrbitalBasis<double, GTOShell>& spin_orbital_basis) {
                return spin_orbital_basis.quantize(Operator::Overlap());
            },
            "Return the overlap operator expressed in this spinor basis.")


        .def(
            "quantizeCoulombRepulsionOperator",
            [](const RSpinOrbitalBasis<double, GTOShell>& spin_orbital_basis) {
                return spin_orbital_basis.quantize(Operator::Coulomb());
            },
            "Return the Coulomb repulsion operator expressed in this spinor basis.")


        /*
         *  MARK: Mulliken partitioning
         */

        .def(
            "mullikenPartitioning",
            [](const RSpinOrbitalBasis<double, GTOShell>& spin_orbital_basis, const std::function<bool(const GTOShell&)>& selector) {
                spin_orbital_basis.mullikenPartitioning(selector);
            },
            py::arg("selector"),
            "A `RMullikenPartitioning` for the AOs selected by the supplied selector function.")


        /*
         *  MARK: Conforming to `BasisTransformable`
         */

        .def(
            "rotate",
            [](RSpinOrbitalBasis<double, GTOShell>& spin_orbital_basis, const RTransformation<double>& U) {
                spin_orbital_basis.rotate(U);
            },
            py::arg("U"),
            "In-place apply the basis rotation.")

        .def(
            "rotated",
            [](const RSpinOrbitalBasis<double, GTOShell>& spin_orbital_basis, const RTransformation<double>& U) {
                spin_orbital_basis.rotated(U);
            },
            py::arg("U"),
            "Apply the basis rotation and return the result.")

        .def(
            "transform",
            [](RSpinOrbitalBasis<double, GTOShell>& spin_orbital_basis, const RTransformation<double>& T) {
                spin_orbital_basis.transform(T);
            },
            py::arg("T"),
            "In-place apply the basis transformation.")

        .def(
            "transformed",
            [](const RSpinOrbitalBasis<double, GTOShell>& spin_orbital_basis, const RTransformation<double>& T) {
                spin_orbital_basis.transformed(T);
            },
            py::arg("T"),
            "Apply the basis transformation and return the result.");
}


}  // namespace gqcpy
