// This file is part of GQCG-gqcp.
// 
// Copyright (C) 2017-2018  the GQCG developers
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
#define BOOST_TEST_MODULE "ERNewtonLocalizer"


#include "Localization/ERNewtonLocalizer.hpp"

#include "HamiltonianParameters/HamiltonianParameters.hpp"


#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain


BOOST_AUTO_TEST_CASE ( localization_index_raises ) {

    // Check if the Edmiston-Ruedenberg localization index is raised after a localization procedure
    GQCP::Molecule h2o ("../tests/data/h2o.xyz");
    size_t N_P = h2o.get_N()/2;

    auto mol_ham_par = GQCP::HamiltonianParameters::Molecular(h2o, "STO-3G");  // in AOBasis
    mol_ham_par.LowdinOrthonormalize();  // in the Löwdin basis


    double D_before = mol_ham_par.calculateEdmistonRuedenbergLocalizationIndex(N_P);

    GQCP::ERNewtonLocalizer localizer (N_P, 1.0e-04);
    localizer.localize(mol_ham_par);  // now the Hamiltonian parameters are in the localized basis

    double D_after = mol_ham_par.calculateEdmistonRuedenbergLocalizationIndex(N_P);

    BOOST_CHECK(D_after > D_before);
}
