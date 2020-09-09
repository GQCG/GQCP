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
#include "Basis/SpinorBasis/USpinorBasis.hpp"
#include "Molecule/Molecule.hpp"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>


namespace py = pybind11;


namespace gqcpy {


void bindUSpinorBasis(py::module& module) {
    py::class_<GQCP::USpinorBasis<double, GQCP::GTOShell>>(module, "USpinorBasis", "A class that represents a real, unrestricted spinor basis with underlying GTO shells.")

        // CONSTRUCTORS

        .def(py::init<const GQCP::Molecule&, const std::string&>(),
             py::arg("molecule"),
             py::arg("basisset_name"))

        .def_static(
            "FromRestricted",
            [](const GQCP::RSpinorBasis<double, GQCP::GTOShell>& r_spinor_basis) {
                return GQCP::USpinorBasis<double, GQCP::GTOShell>::FromRestricted(r_spinor_basis);
            },
            py::arg("r_spinor_basis"),
            "Create an unrestricted spinor basis from a restricted spinor basis, leading to alpha- and beta- coefficient matrices that are equal.")


        // PUBLIC METHODS

        .def(
            "calculateAtomicSpinZ",
            [](const GQCP::USpinorBasis<double, GQCP::GTOShell>& spinor_basis, const std::vector<size_t>& ao_list, const GQCP::Spin sigma) {
                return spinor_basis.calculateAtomicSpinZ(ao_list, sigma);
            },
            py::arg("ao_list"),
            py::arg("sigma"),
            "Return the atomic spin-z operator for a given set of AO indices.")

        .def(
            "calculateMullikenOperator",
            [](GQCP::USpinorBasis<double, GQCP::GTOShell>& spinor_basis, const std::vector<size_t>& ao_list, const GQCP::Spin sigma) {
                return spinor_basis.calculateMullikenOperator(ao_list, sigma);
            },
            py::arg("ao_list"),
            py::arg("sigma"),
            "Return the Mulliken operator for a set of given AO indices.")

        .def(
            "coefficientMatrix",
            [](const GQCP::USpinorBasis<double, GQCP::GTOShell>& spinor_basis, const GQCP::Spin sigma) {
                return spinor_basis.coefficientMatrix(sigma);
            },
            py::arg("sigma"),
            "Return the coefficient matrix for the requested component, i.e. the matrix of the expansion coefficients of the requested component of the spinors in terms of its underlying scalar basis")

        .def(
            "numberOfSpinors",
            [](const GQCP::USpinorBasis<double, GQCP::GTOShell>& spinor_basis) {
                return spinor_basis.numberOfSpinors();
            },
            "Return the total number of spinors/spin-orbitals that this spinor basis describes")

        .def(
            "numberOfSpinors",
            [](const GQCP::USpinorBasis<double, GQCP::GTOShell>& spinor_basis, const GQCP::Spin sigma) {
                return spinor_basis.numberOfSpinors(sigma);
            },
            py::arg("sigma"),
            "Return the total number of sigma-spinors/spin-orbitals that this spinor basis describes")

        .def(
            "isOrthornomal",
            [](const GQCP::USpinorBasis<double, GQCP::GTOShell>& spinor_basis, const GQCP::Spin sigma, const double precision) {
                return spinor_basis.isOrthonormal(sigma, precision);
            },
            py::arg("sigma"),
            py::arg("precision") = 1.0e-08,
            "Return if this spinor basis for the requested component is orthonormal within the given precision")

        .def(
            "isOrthornomal",
            [](const GQCP::USpinorBasis<double, GQCP::GTOShell>& spinor_basis, const double precision) {
                return spinor_basis.isOrthonormal(precision);
            },
            py::arg("precision") = 1.0e-08,
            "Return if this spinor basis is orthonormal within the given precision")

        .def(
            "lowdinOrthonormalize",
            [](GQCP::USpinorBasis<double, GQCP::GTOShell>& spinor_basis) {
                spinor_basis.lowdinOrthonormalize();
            },
            "Transform the spinor basis to the 'Löwdin basis', which is the orthonormal basis that we transform to with T = S^{-1/2}, where S is the current overlap matrix.")

        .def(
            "overlap",
            [](const GQCP::USpinorBasis<double, GQCP::GTOShell>& spinor_basis, const GQCP::Spin sigma) {
                return spinor_basis.overlap(sigma);
            },
            py::arg("sigma"),
            "Return the overlap (one-electron) operator of the requested component of this spinor basis")

        .def(
            "quantizeCoulombRepulsionOperator",
            [](const GQCP::USpinorBasis<double, GQCP::GTOShell>& spinor_basis, const GQCP::Spin sigma) {
                return spinor_basis.quantize(GQCP::Operator::Coulomb(), sigma);
            },
            py::arg("sigma"),
            "Return the Coulomb repulsion operator expressed in this spinor basis.")

        .def(
            "quantizeDipoleOperator",
            [](const GQCP::USpinorBasis<double, GQCP::GTOShell>& spinor_basis, const GQCP::Vector<double, 3>& origin, const GQCP::Spin sigma) {
                return spinor_basis.quantize(GQCP::Operator::ElectronicDipole(origin), sigma);
            },
            py::arg("spin"),
            py::arg("origin") = GQCP::Vector<double, 3>::Zero(),
            "Return the electronic dipole operator expressed in this spinor basis.")

        .def(
            "quantizeKineticOperator",
            [](const GQCP::USpinorBasis<double, GQCP::GTOShell>& spinor_basis, const GQCP::Spin sigma) {
                return spinor_basis.quantize(GQCP::Operator::Kinetic(), sigma);
            },
            py::arg("sigma"),
            "Return the kinetic energy operator expressed in this spinor basis.")

        .def(
            "quantizeNuclearAttractionOperator",
            [](const GQCP::USpinorBasis<double, GQCP::GTOShell>& spinor_basis, const GQCP::Molecule& molecule, const GQCP::Spin sigma) {
                return spinor_basis.quantize(GQCP::Operator::NuclearAttraction(molecule), sigma);
            },
            py::arg("molecule"),
            py::arg("sigma"),
            "Return the nuclear attraction operator expressed in this spinor basis.")

        .def(
            "quantizeOverlapOperator",
            [](const GQCP::USpinorBasis<double, GQCP::GTOShell>& spinor_basis, const GQCP::Spin sigma) {
                return spinor_basis.quantize(GQCP::Operator::Overlap(), sigma);
            },
            py::arg("sigma"),
            "Return the overlap operator expressed in this spinor basis.")

        .def(
            "rotate",
            [](GQCP::USpinorBasis<double, GQCP::GTOShell>& spinor_basis, const Eigen::MatrixXd& U, const GQCP::Spin sigma) {
                spinor_basis.rotate(GQCP::TransformationMatrix<double>(U), sigma);
            },
            py::arg("U"),
            py::arg("sigma"),
            "Rotate the spinor basis of the requested component to another one using the given unitary transformation matrix.")

        .def(
            "transform",
            [](GQCP::USpinorBasis<double, GQCP::GTOShell>& spinor_basis, const Eigen::MatrixXd& T, const GQCP::Spin sigma) {
                spinor_basis.transform(GQCP::TransformationMatrix<double>(T), sigma);
            },
            py::arg("T"),
            py::arg("sigma"),
            "Transform the spinor basis for one component to another one using the given transformation matrix.")

        .def(
            "transform",
            [](GQCP::USpinorBasis<double, GQCP::GTOShell>& spinor_basis, const Eigen::MatrixXd& T_alpha, const Eigen::MatrixXd& T_beta) {
                spinor_basis.transform(GQCP::SpinResolvedTransformationMatrix<double>(T_alpha, T_beta));
            },
            py::arg("T_alpha"),
            py::arg("T_beta"),
            "Transform the spinor basis to another one using the given transformation matrices.");
}


}  // namespace gqcpy
