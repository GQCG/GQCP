#define BOOST_TEST_MODULE "HamiltonianParameters_constructors"


#include "HamiltonianParameters/HamiltonianParameters_constructors.hpp"


#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain



BOOST_AUTO_TEST_CASE ( constructMolecularHamiltonianParameters ) {

    // Set up a basis
    GQCG::Molecule water ("../tests/data/h2o.xyz");
    auto ao_basis_sptr = std::make_shared<GQCG::AOBasis>(water, "STO-3G");

    
    // Check if we can construct the molecular Hamiltonian parameters
    auto mol_ham_par = GQCG::constructMolecularHamiltonianParameters(ao_basis_sptr);
}
