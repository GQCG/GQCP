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

#define BOOST_TEST_MODULE "ERNewtonLocalizer"

#include <boost/test/unit_test.hpp>

#include "Mathematical/Optimization/Minimization/IterativeIdentitiesHessianModifier.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "QCMethod/OrbitalOptimization/Localization/ERNewtonLocalizer.hpp"


/**
 *  Check if the Edmiston-Ruedenberg localization index is raised after a localization procedure.
 * 
 *  The test system is H2O in an STO-3G basisset.
 */
BOOST_AUTO_TEST_CASE(localization_index_raises) {

    // Prepare the molecular Hamiltonian in the Löwdin-basis.
    const auto molecule = GQCP::Molecule::ReadXYZ("data/h2o.xyz");
    const auto N_P = molecule.numberOfElectronPairs();

    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis {molecule, "STO-3G"};
    spinor_basis.lowdinOrthonormalize();
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, molecule);  // in the Löwdin basis


    // Do an Edmiston-Ruedenberg localization and keep track of the value of the ER-index before and after.
    const auto orbital_space = GQCP::OrbitalSpace::Implicit({{GQCP::OccupationType::k_occupied, N_P}});  // N_P occupied spatial orbitals
    double D_before = sq_hamiltonian.calculateEdmistonRuedenbergLocalizationIndex(orbital_space);

    auto hessian_modifier = std::make_shared<GQCP::IterativeIdentitiesHessianModifier>();
    GQCP::ERNewtonLocalizer localizer {orbital_space, hessian_modifier, 1.0e-04};
    localizer.optimize(spinor_basis, sq_hamiltonian);  // if converged, the Hamiltonian is in the localized basis

    double D_after = sq_hamiltonian.calculateEdmistonRuedenbergLocalizationIndex(orbital_space);

    BOOST_CHECK(D_after > D_before);
}
