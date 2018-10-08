#define BOOST_TEST_MODULE "AP1roG"


#include "AP1roG/AP1roGPSESolver.hpp"
#include "HamiltonianParameters/HamiltonianParameters_constructors.hpp"


#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain


BOOST_AUTO_TEST_CASE ( constructor ) {

    // Test a correct constructor
    GQCG::Molecule h2 ("../tests/data/h2_szabo.xyz");
    auto ao_basis = std::make_shared<GQCG::AOBasis>(h2, "STO-3G");
    auto mol_ham_par = GQCG::constructMolecularHamiltonianParameters(ao_basis);
    GQCG::AP1roGPSESolver ap1rog_pse_solver (h2, mol_ham_par);


    // Test a faulty constructor
    GQCG::Molecule h2_ion ("../tests/data/h2_szabo.xyz", +1);
    BOOST_CHECK_THROW(GQCG::AP1roGPSESolver(h2_ion, mol_ham_par), std::invalid_argument);  // we can use the same Hamiltonian parameters for molecule and ion
}
