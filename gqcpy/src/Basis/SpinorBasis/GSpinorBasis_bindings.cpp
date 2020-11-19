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

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>


namespace py = pybind11;


namespace gqcpy {


void bindGSpinorBasis(py::module& module) {
    py::class_<GQCP::GSpinorBasis<double, GQCP::GTOShell>>(module, "GSpinorBasis", "A class that represents a real, (generalized) spinor basis with underlying GTO shells.")

        // CONSTRUCTORS

        .def(py::init<const GQCP::Molecule&, const std::string&>(),
             py::arg("molecule"),
             py::arg("basisset_name"))


        .def_static("FromRestricted",
                    [](const GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell>& r_spinor_basis) {
                        return GQCP::GSpinorBasis<double, GQCP::GTOShell>::FromRestricted(r_spinor_basis);
                    })

        // INHERITED METHODS

        .def(
            "expansion",
            [](const GQCP::GSpinorBasis<double, GQCP::GTOShell>& spinor_basis) {
                return spinor_basis.expansion();
            },
            "Return the transformation matrix between the scalar bases and the current spinors.")

        .def("expansion",
             [](const GQCP::GSpinorBasis<double, GQCP::GTOShell>& spinor_basis, const GQCP::Spin sigma) {
                 return spinor_basis.expansion(sigma);
             })

        .def(
            "isOrthonormal",
            [](const GQCP::GSpinorBasis<double, GQCP::GTOShell>& spinor_basis, const double precision) {
                return spinor_basis.isOrthonormal(precision);
            },
            py::arg("precision") = 1.0e-08,
            "Return if this spinor basis is orthonormal within the given precision")

        .def(
            "lowdinOrthonormalization",
            [](const GQCP::GSpinorBasis<double, GQCP::GTOShell>& spinor_basis) {
                return spinor_basis.lowdinOrthonormalization();
            },
            "Return the transformation matrix to the Löwdin basis: T = S_current^{-1/2}")

        .def(
            "lowdinOrthonormalize",
            [](GQCP::GSpinorBasis<double, GQCP::GTOShell>& spinor_basis) {
                spinor_basis.lowdinOrthonormalize();
            },
            "Transform the spinor basis to the 'Löwdin basis', which is the orthonormal basis that we transform to with T = S^{-1/2}, where S is the current overlap matrix.")

        .def(
            "overlap",
            [](const GQCP::GSpinorBasis<double, GQCP::GTOShell>& spinor_basis) {
                return spinor_basis.overlap();
            },
            "Return the overlap (one-electron) operator of this restricted spinor basis")

        .def(
            "rotate",
            [](GQCP::GSpinorBasis<double, GQCP::GTOShell>& spinor_basis, const Eigen::MatrixXd& U) {
                spinor_basis.rotate(GQCP::GTransformation<double>(U));
            },
            py::arg("U"),
            "Rotate the spinor basis to another one using the given unitary transformation matrix.")

        .def(
            "transform", [](GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell>& spinor_basis, const Eigen::MatrixXd& T) {
                spinor_basis.transform(GQCP::GTransformation<double>(T));
            },
            py::arg("T"), "Transform the current spinor basis using a given transformation matrix")


        // PUBLIC METHODS

        .def(
            "expansion", [](const GQCP::GSpinorBasis<double, GQCP::GTOShell>& spinor_basis, const GQCP::Spin sigma) {
                return spinor_basis.expansion(sigma);
            },
            py::arg("sigma"), "Return the coefficient matrix for the requested spin component, i.e. the matrix of the expansion coefficients of the requested components of the spinors in terms of its underlying scalar basis")

        .def(
            "numberOfCoefficients", [](const GQCP::GSpinorBasis<double, GQCP::GTOShell>& spinor_basis, const GQCP::Spin sigma) {
                return spinor_basis.numberOfCoefficients(sigma);
            },
            py::arg("sigma"), "Return the number of coefficients that are used for the expansion of the requested spin-component of a spinor")

        .def(
            "numberOfSpinors", [](const GQCP::GSpinorBasis<double, GQCP::GTOShell>& spinor_basis) {
                return spinor_basis.numberOfSpinors();
            },
            "Return the number of spinors that 'are' in this generalized spinor basis.")

        .def(
            "quantizeCoulombRepulsionOperator", [](const GQCP::GSpinorBasis<double, GQCP::GTOShell>& spinor_basis) {
                return spinor_basis.quantize(GQCP::Operator::Coulomb());
            },
            "Return the Coulomb operator expressed in this spinor basis.")

        .def(
            "quantizeKineticOperator", [](const GQCP::GSpinorBasis<double, GQCP::GTOShell>& spinor_basis) {
                return spinor_basis.quantize(GQCP::Operator::Kinetic());
            },
            "Return the kinetic energy operator expressed in this spinor basis.")

        .def(
            "quantizeNuclearAttractionOperator", [](const GQCP::GSpinorBasis<double, GQCP::GTOShell>& spinor_basis, const GQCP::Molecule& molecule) {
                return spinor_basis.quantize(GQCP::Operator::NuclearAttraction(molecule));
            },
            "Return the nuclear attraction operator expressed in this spinor basis.")

        .def(
            "quantizeOverlapOperator", [](const GQCP::GSpinorBasis<double, GQCP::GTOShell>& spinor_basis) {
                return spinor_basis.quantize(GQCP::Operator::Overlap());
            },
            "Return the overlap operator expressed in this spinor basis.")

        .def(
            "quantizeSpinOperator", [](const GQCP::GSpinorBasis<double, GQCP::GTOShell>& spinor_basis) {
                return spinor_basis.quantize(GQCP::Operator::ElectronicSpin());
            },
            "Return the electronic spin operator expressed in this spinor basis.")

        .def(
            "transform", [](GQCP::GSpinorBasis<double, GQCP::GTOShell>& spinor_basis, const Eigen::MatrixXd& T_matrix) {
                const GQCP::GTransformation<double> T(T_matrix);
                spinor_basis.transform(T);
            },
            "Transform the current spinor basis using a given transformation matrix.");
}


}  // namespace gqcpy
