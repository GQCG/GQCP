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

#define BOOST_TEST_MODULE "QCMethod_GHF_stability_test"

#include <boost/test/unit_test.hpp>

#include "Basis/SpinorBasis/GSpinorBasis.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "QCMethod/HF/GHF/GHF.hpp"
#include "QCMethod/HF/GHF/GHFSCFSolver.hpp"
#include "QCMethod/HF/GHF/GHFStabilityChecks.hpp"
#include "QCModel/HF/Stability/StabilityMatrices/GHFStabilityMatrix.hpp"

/**
 *  Check if the plain GHF SCF solver finds a correct solution.
 * 
 *  The system of interest is a H3-triangle, 1 bohr apart and the reference implementation was done by @xdvriend.
 */
BOOST_AUTO_TEST_CASE(H3_stability_test_1) {

    // Set up a general spinor basis to obtain a spin-blocked second-quantized molecular Hamiltonian.
    const auto molecule = GQCP::Molecule::HRingFromDistance(3, 1.0);  // H3-triangle, 1 bohr apart
    //const auto molecule = GQCP::Molecule::HChain(2, 1.0);
    const auto N = molecule.numberOfElectrons();

    const GQCP::GSpinorBasis<double, GQCP::GTOShell> g_spinor_basis {molecule, "STO-3G"};
    const auto S = g_spinor_basis.overlap().parameters();

    const auto sq_hamiltonian = GQCP::GSQHamiltonian<double>::Molecular(g_spinor_basis, molecule);

    auto environment = GQCP::GHFSCFEnvironment<double>::WithCoreGuess(N, sq_hamiltonian, S);

    auto solver = GQCP::GHFSCFSolver<double>::Plain(1.0e-08, 3000);
    const auto qc_structure = GQCP::QCMethod::GHF<double>().optimize(solver, environment);
    auto ghf_parameters = qc_structure.groundStateParameters();

    ghf_parameters.stabilityProperties().print();

    GQCP::QCMethod::GHFStabilityChecks<double>().internalStabilityCheck(ghf_parameters, sq_hamiltonian);

    ghf_parameters.stabilityProperties().print();

    GQCP::QCMethod::GHFStabilityChecks<double>().externalStabilityCheck(ghf_parameters, sq_hamiltonian);

    ghf_parameters.stabilityProperties().print();
}