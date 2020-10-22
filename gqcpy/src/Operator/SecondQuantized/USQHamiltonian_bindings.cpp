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
#include "Operator/SecondQuantized/USQHamiltonian.hpp"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>


namespace py = pybind11;


namespace gqcpy {


void bindUSQHamiltonian(py::module& module) {
    py::class_<GQCP::USQHamiltonian<double>>(module, "USQHamiltonian", "A class that represents a real, second-quantized unrestricted Hamiltonian.")

        // CONSTRUCTORS

        .def_static(
            "Molecular",
            [](const GQCP::USpinorBasis<double, GQCP::GTOShell>& u_spinor_basis, const GQCP::Molecule& molecule) {
                return GQCP::USQHamiltonian<double>::Molecular(u_spinor_basis, molecule);
            },
            py::arg("u_spinor_basis"),
            py::arg("molecule"),
            "Construct the unrestricted molecular Hamiltonian in a given unrestricted spin-orbital basis.")

        // PUBLIC METHODS

        .def(
            "areSpinHamiltoniansOfSameDimension",
            [](const GQCP::USQHamiltonian<double>& usq_hamiltonian) {
                return usq_hamiltonian.areSpinHamiltoniansOfSameDimension();
            },
            "Return whether the spin components of the unrestricted Hamiltonian have the same dimension.")

        .def(
            "constrain",
            [](GQCP::USQHamiltonian<double>& usq_hamiltonian, const GQCP::ScalarRSQOneElectronOperator<double>& one_electron_op, const double lambda, const GQCP::Spin sigma) {
                return usq_hamiltonian.constrain(one_electron_op, lambda, sigma);
            },
            py::arg("one_electron_op"),
            py::arg("lambda"),
            py::arg("sigma"),
            "Constrain spin component sigma of the unrestricted Hamiltonian with a certain one electron operator and Lagrangian lambda.")

        .def(
            "numberOfOrbitals",
            [](const GQCP::USQHamiltonian<double>& usq_hamiltonian) {
                return usq_hamiltonian.numberOfOrbitals();
            },
            "Return the total number of orbitals in the USQ Hamiltonian.")

        .def(
            "rotate",
            [](GQCP::USQHamiltonian<double>& usq_hamiltonian, const Eigen::MatrixXd& U) {
                usq_hamiltonian.rotate(GQCP::TransformationMatrix<double> {U});
            },
            py::arg("U"),
            "Rotate both components of the USQ Hamiltonian with matrix U.")

        .def(
            "rotate",
            [](GQCP::USQHamiltonian<double>& usq_hamiltonian, const Eigen::MatrixXd& U, const GQCP::Spin sigma) {
                usq_hamiltonian.rotate(GQCP::TransformationMatrix<double> {U}, sigma);
            },
            py::arg("U"),
            py::arg("sigma"),
            "Rotate the spin sigma component of the USQ Hamiltonian with matrix U.")

        .def(
            "rotate",
            [](GQCP::USQHamiltonian<double>& usq_hamiltonian, const GQCP::JacobiRotationParameters& jacobi_rotation_parameters) {
                usq_hamiltonian.rotate(jacobi_rotation_parameters);
            },
            py::arg("jacobi_rotation_parameters"),
            "Rotate both components of the USQ Hamiltonian using jacobi parameters.")

        .def(
            "rotate",
            [](GQCP::USQHamiltonian<double>& usq_hamiltonian, const GQCP::JacobiRotationParameters& jacobi_rotation_parameters, const GQCP::Spin sigma) {
                usq_hamiltonian.rotate(jacobi_rotation_parameters);
            },
            py::arg("jacobi_rotation_parameters"),
            py::arg("sigma"),
            "Rotate the spin sigma component of the USQ Hamiltonian using jacobi parameters.")

        .def(
            "spinHamiltonian",
            [](const GQCP::USQHamiltonian<double>& usq_hamiltonian, const GQCP::Spin sigma) {
                return usq_hamiltonian.spinHamiltonian(sigma);
            },
            py::arg("sigma"),
            "Return the spin sigma component of the USQ Hamiltonian")

        .def(
            "transform",
            [](GQCP::USQHamiltonian<double>& usq_hamiltonian, const Eigen::MatrixXd& T) {
                usq_hamiltonian.transform(GQCP::TransformationMatrix<double> {T});
            },
            py::arg("U"),
            "Transform both components of the USQ Hamiltonian with matrix U.")

        .def(
            "transform",
            [](GQCP::USQHamiltonian<double>& usq_hamiltonian, const Eigen::MatrixXd& T, const GQCP::Spin sigma) {
                usq_hamiltonian.transform(GQCP::TransformationMatrix<double> {T}, sigma);
            },
            py::arg("U"),
            py::arg("sigma"),
            "Transform the spin sigma component of the USQ Hamiltonian with matrix U.")

        .def(
            "twoElectronContributionsMixed",
            [](const GQCP::USQHamiltonian<double>& usq_hamiltonian) {
                return usq_hamiltonian.twoElectronContributionsMixed();
            },
            "Return the contributions to the mixed alpha & beta two-electron part of the unrestricted Hamiltonian.")

        .def(
            "twoElectronMixed",
            [](const GQCP::USQHamiltonian<double>& usq_hamiltonian) {
                return usq_hamiltonian.twoElectronMixed();
            },
            "Return the total contributions to the mixed alpha & beta two-electron part of the unrestricted Hamiltonian.");
}


}  // namespace gqcpy
