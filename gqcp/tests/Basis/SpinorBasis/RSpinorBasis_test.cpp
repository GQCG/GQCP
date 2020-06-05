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

#define BOOST_TEST_MODULE "RSpinorBasis_test"

#include <boost/test/unit_test.hpp>

#include "Basis/SpinorBasis/RSpinorBasis.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "QCMethod/HF/RHF/DiagonalRHFFockMatrixObjective.hpp"
#include "QCMethod/HF/RHF/RHF.hpp"
#include "QCMethod/HF/RHF/RHFSCFSolver.hpp"


/**
 *  Check if the orbitals in an AO basis are not orthonormal, but after a transformation to the canonical RHF orbitals, they are.
 */
BOOST_AUTO_TEST_CASE(RHF_orbitals_are_orthonormal) {

    // The orbitals in an AO basis are not orthonormal.
    const auto molecule = GQCP::Molecule::ReadXYZ("data/h2o.xyz");
    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis {molecule, "STO-3G"};
    BOOST_CHECK(!spinor_basis.isOrthonormal());


    // The orbitals in the RHF basis should be orthonormal.
    const auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, molecule);  // in the scalar/AO basis

    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(molecule.numberOfElectrons(), sq_hamiltonian, spinor_basis.overlap().parameters());
    auto plain_rhf_scf_solver = GQCP::RHFSCFSolver<double>::Plain();
    const GQCP::DiagonalRHFFockMatrixObjective<double> objective {sq_hamiltonian};
    const auto rhf_parameters = GQCP::QCMethod::RHF<double>().optimize(objective, plain_rhf_scf_solver, rhf_environment).groundStateParameters();

    spinor_basis.transform(rhf_parameters.coefficientMatrix());
    BOOST_CHECK(spinor_basis.isOrthonormal());
}


#include "Mathematical/Grid/CubicGrid.hpp"


/**
 *  Check if the orbitals in an AO basis are not orthonormal, but after a transformation to the canonical RHF orbitals, they are
 */
BOOST_AUTO_TEST_CASE(sandbox) {

    // The orbitals in an AO basis are not orthonormal
    const auto molecule = GQCP::Molecule::ReadXYZ("data/h2o.xyz");
    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis {molecule, "STO-3G"};


    // The orbitals in the RHF basis should be orthonormal
    const auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, molecule);  // in the scalar/AO basis

    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(molecule.numberOfElectrons(), sq_hamiltonian, spinor_basis.overlap().parameters());
    auto plain_rhf_scf_solver = GQCP::RHFSCFSolver<double>::Plain();
    const GQCP::DiagonalRHFFockMatrixObjective<double> objective {sq_hamiltonian};
    const auto rhf_parameters = GQCP::QCMethod::RHF<double>().optimize(objective, plain_rhf_scf_solver, rhf_environment).groundStateParameters();

    spinor_basis.transform(rhf_parameters.coefficientMatrix());


    const auto spin_orbitals = spinor_basis.spinOrbitals();
    // const auto lumo_index = 2 * rhf_parameters.lumoIndex();  // convert spatial orbital to spin-orbital index
    // const auto homo_index = 2 * rhf_parameters.homoIndex();
    const size_t homo_alpha_index = 8;
    const size_t lumo_alpha_index = 10;


    const auto& alpha_homo = spin_orbitals[homo_alpha_index].component(GQCP::Spin::alpha);
    const auto& alpha_lumo = spin_orbitals[lumo_alpha_index].component(GQCP::Spin::alpha);


    const auto grid = GQCP::CubicGrid::Centered(molecule.nuclearFramework().nucleiAsVector()[0].position(), 50, 0.1);
    const auto lumo_field = grid.evaluate(alpha_homo);
    const auto homo_field = grid.evaluate(alpha_lumo);


    grid.writeToCubeFile(homo_field, "homo.cube", molecule);
    grid.writeToCubeFile(lumo_field, "lumo.cube", molecule);
}
