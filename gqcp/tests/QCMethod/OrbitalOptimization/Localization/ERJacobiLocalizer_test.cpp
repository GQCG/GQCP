// This file is part of GQCG-gqcp.
// 
// Copyright (C) 2017-2019  the GQCG developers
// 
// GQCG-gqcp is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// GQCG-gqcp is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-gqcp.  If not, see <http://www.gnu.org/licenses/>.
// 
#define BOOST_TEST_MODULE "ERJacobiLocalizer"

#include <boost/test/unit_test.hpp>

#include "QCMethod/OrbitalOptimization/Localization/ERJacobiLocalizer.hpp"


BOOST_AUTO_TEST_CASE ( localization_index_raises ) {

    // Check if the Edmiston-Ruedenberg localization index is raised after a localization procedure
    auto h2o = GQCP::Molecule::ReadXYZ("data/h2o.xyz");
    size_t N_P = h2o.numberOfElectrons()/2;

    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis (h2o, "STO-3G");
    spinor_basis.lowdinOrthonormalize();
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, h2o);  // in the Löwdin basis


    double D_before = sq_hamiltonian.calculateEdmistonRuedenbergLocalizationIndex(N_P);

    GQCP::ERJacobiLocalizer localizer (N_P, 1.0e-04);
    localizer.optimize(spinor_basis, sq_hamiltonian);  // now the Hamiltonian is in the localized basis

    double D_after = sq_hamiltonian.calculateEdmistonRuedenbergLocalizationIndex(N_P);

    BOOST_CHECK(D_after > D_before);
}
