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

#define BOOST_TEST_MODULE "expectation_values"

#include <boost/test/unit_test.hpp>

#include "Basis/Integrals/Interfaces/LibintInterfacer.hpp"
#include "Basis/transform.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "Operator/SecondQuantized/USQHamiltonian.hpp"
#include "Processing/Properties/expectation_values.hpp"
#include "QCMethod/CI/HamiltonianBuilder/DOCI.hpp"
#include "QCMethod/CI/HamiltonianBuilder/FCI.hpp"
#include "QCMethod/HF/RHF/DiagonalRHFFockMatrixObjective.hpp"
#include "QCMethod/HF/RHF/RHF.hpp"
#include "QCMethod/HF/RHF/RHFSCFSolver.hpp"
#include "Utilities/units.hpp"


/**
 *  Check that the total RHF Mulliken population of N2 is 14.
 */
BOOST_AUTO_TEST_CASE(mulliken_N2_STO_3G) {

    // Initialize the molecular Hamiltonian for N2 in the Löwdin-orthonormalized basis.
    const GQCP::Nucleus N1 {7, 0.0, 0.0, 0.0};
    const GQCP::Nucleus N2 {7, 0.0, 0.0, GQCP::units::angstrom_to_bohr(1.134)};  // from CCCBDB, STO-3G geometry
    const GQCP::Molecule molecule {{N1, N2}};
    const auto N = molecule.numberOfElectrons();

    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis {molecule, "STO-3G"};
    spinor_basis.lowdinOrthonormalize();
    const auto K = spinor_basis.numberOfSpatialOrbitals();

    const auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, molecule);


    // Calculate the Mulliken population operator in this spinor basis.

    // To calculate the total Mulliken population operator, we have to include all basis functions.
    std::vector<size_t> basis_functions;
    basis_functions.reserve(K);
    for (size_t i = 0; i < K; i++) {
        basis_functions.push_back(i);
    }
    const auto mulliken_op = spinor_basis.calculateMullikenOperator(basis_functions);  // 'op' for 'operator'


    // Create the RHF 1-DM for N2 and check the total Mulliken operator.
    const auto D = GQCP::QCModel::RHF<double>::calculateOrthonormalBasis1DM(K, N);

    const auto mulliken_population = mulliken_op.calculateExpectationValue(D)(0);
    BOOST_CHECK(std::abs(mulliken_population - N) < 1.0e-12);
}
