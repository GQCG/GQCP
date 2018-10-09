#define BOOST_TEST_MODULE "DOCI"


#include "HamiltonianBuilder/DOCI.hpp"
#include "HamiltonianParameters/HamiltonianParameters_constructors.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain



BOOST_AUTO_TEST_CASE ( DOCI_constructor  ) {
    // Create an AOBasis
    GQCG::Molecule water ("../tests/data/h2o.xyz");
    auto ao_basis_ptr = std::make_shared<GQCG::AOBasis>(water, "STO-3G");


    // Create random HamiltonianParameters from One- and TwoElectronOperators (and a transformation matrix) with compatible dimensions
    size_t K = ao_basis_ptr->get_number_of_basis_functions();
    GQCG::OneElectronOperator S (Eigen::MatrixXd::Random(K, K));
    GQCG::OneElectronOperator H_core (Eigen::MatrixXd::Random(K, K));
    Eigen::Tensor<double, 4> g_tensor (K, K, K, K);
    g_tensor.setRandom();
    GQCG::TwoElectronOperator g (g_tensor);
    Eigen::MatrixXd C = Eigen::MatrixXd::Random(K, K);
    GQCG::HamiltonianParameters random_hamiltonian_parameters (ao_basis_ptr, S, H_core, g, C);

    // Create a compatible Fock space
    GQCG::FockSpace fock_space (K, 3);

    // Check if a correct constructor works and public methods don't fail
    BOOST_CHECK_NO_THROW(GQCG::DOCI random_doci(random_hamiltonian_parameters, fock_space));

    // Check if faulty constructor parameters throw an error
    // Create an incompatible Fock space
    GQCG::FockSpace fock_space_2 (K+1, 3);
    BOOST_CHECK_THROW(GQCG::DOCI random_doci_2(random_hamiltonian_parameters,fock_space_2), std::invalid_argument);
}


BOOST_AUTO_TEST_CASE ( DOCI_public_methods ) {
    // Create an AOBasis
    GQCG::Molecule water ("../tests/data/h2o.xyz");
    auto ao_basis_ptr = std::make_shared<GQCG::AOBasis>(water, "STO-3G");


    // Create random HamiltonianParameters from One- and TwoElectronOperators (and a transformation matrix) with compatible dimensions
    size_t K = ao_basis_ptr->get_number_of_basis_functions();
    GQCG::OneElectronOperator S (Eigen::MatrixXd::Random(K, K));
    GQCG::OneElectronOperator H_core (Eigen::MatrixXd::Random(K, K));
    Eigen::Tensor<double, 4> g_tensor (K, K, K, K);
    g_tensor.setRandom();
    GQCG::TwoElectronOperator g (g_tensor);
    Eigen::MatrixXd C = Eigen::MatrixXd::Random(K, K);
    GQCG::HamiltonianParameters random_hamiltonian_parameters (ao_basis_ptr, S, H_core, g, C);

    // Create a compatible Fock space
    GQCG::FockSpace fock_space (K, 3);

    // Create DOCI module
    GQCG::DOCI random_doci (random_hamiltonian_parameters,fock_space);

    // Test the public DOCI methods
    Eigen::VectorXd x = random_doci.calculateDiagonal();
    BOOST_CHECK_NO_THROW(random_doci.constructHamiltonian());
    BOOST_CHECK_NO_THROW(random_doci.matrixVectorProduct(x));
}


BOOST_AUTO_TEST_CASE ( DOCI_beh_cation_klaas_dense ) {

   // Klaas' reference DOCI energy for BeH+ (obtained through Caitlin)
   double reference_doci_energy = -14.8782216937;

   // Do a DOCI calculation based on a given FCIDUMP file
   // Create the Hamiltonian Parameters
   auto ham_par = GQCG::readFCIDUMPFile("../tests/data/beh_cation_631g_caitlin.FCIDUMP");

   // The species contains 4 electrons and 16 basis functions, this requires a single Fock Space of 16 orbitals and 2 electrons
   GQCG::FockSpace fock_space (16, 2);  // dim = 120

   // Create the DOCI module
   GQCG::DOCI doci (ham_par, fock_space);

   // Retrieve the Hamiltonian matrix
   Eigen::MatrixXd ham = doci.constructHamiltonian();

   // Retrieve the eigenvalues
   Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen_solver (ham);
   auto doci_eigenvalue = eigen_solver.eigenvalues()(0);

   // Calculate the total energy
   double internuclear_repulsion_energy = 1.5900757460937498e+00;  // this comes straight out of the FCIDUMP file
   double test_doci_energy = doci_eigenvalue + internuclear_repulsion_energy;

   BOOST_CHECK(std::abs(test_doci_energy - (reference_doci_energy)) < 1.0e-9);
}


BOOST_AUTO_TEST_CASE ( DOCI_lih_klaas_dense ) {

    // Klaas' reference DOCI energy for LiH (obtained through Caitlin)
    double reference_doci_energy = -8.0029560313;

    // Do a DOCI calculation based on a given FCIDUMP file
    // Create the Hamiltonian Parameters
    auto ham_par = GQCG::readFCIDUMPFile("../tests/data/lih_631g_caitlin.FCIDUMP");

    // The species contains 4 electrons and 16 basis functions, this requires a single Fock Space of 16 orbitals and 2 electrons
    GQCG::FockSpace fock_space (16, 2);  // dim = 120

    // Create the DOCI module
    GQCG::DOCI doci (ham_par,fock_space);

    // Retrieve the Hamiltonian matrix
    Eigen::MatrixXd ham = doci.constructHamiltonian();

    // Retrieve the eigenvalues
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen_solver (ham);
    auto doci_eigenvalue = eigen_solver.eigenvalues()(0);

    // Calculate the total energy
    double internuclear_repulsion_energy = 9.6074293445896852e-01;  // this comes straight out of the FCIDUMP file
    double test_doci_energy = doci_eigenvalue + internuclear_repulsion_energy;

    BOOST_CHECK(std::abs(test_doci_energy - (reference_doci_energy)) < 1.0e-9);
}


BOOST_AUTO_TEST_CASE ( DOCI_li2_klaas_dense ) {

    // Klaas' reference DOCI energy for Li2
    double reference_doci_energy = -15.1153976060;

    // Do a DOCI calculation based on a given FCIDUMP file
    // Create the Hamiltonian Parameters
    auto ham_par = GQCG::readFCIDUMPFile("../tests/data/li2_321g_klaas.FCIDUMP");

    // The species contains 4 electrons and 16 basis functions, this requires a single Fock Space of 16 orbitals and 2 electrons
    GQCG::FockSpace fock_space (18, 3);  // dim = 816

    // Create the DOCI module
    GQCG::DOCI doci (ham_par,fock_space);

    // Retrieve the Hamiltonian matrix
    Eigen::MatrixXd ham = doci.constructHamiltonian();

    // Retrieve the eigenvalues
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen_solver (ham);
    auto doci_eigenvalue = eigen_solver.eigenvalues()(0);

    // Calculate the total energy
    double internuclear_repulsion_energy =  3.0036546888874875e+00;  // this comes straight out of the FCIDUMP file
    double test_doci_energy = doci_eigenvalue + internuclear_repulsion_energy;

    BOOST_CHECK(std::abs(test_doci_energy - (reference_doci_energy)) < 1.0e-9);
}
