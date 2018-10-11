#define BOOST_TEST_MODULE "DOCI"


#include "HamiltonianBuilder/DOCI.hpp"
#include "HamiltonianParameters/HamiltonianParameters_constructors.hpp"
#include "RHF/PlainRHFSCFSolver.hpp"
#include <numopt.hpp>

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain



BOOST_AUTO_TEST_CASE ( DOCI_constructor ) {
    // Create a compatible Fock space
    GQCG::FockSpace fock_space (15, 3);

    // Check if a correct constructor works
    BOOST_CHECK_NO_THROW(GQCG::DOCI doci (fock_space));
}


BOOST_AUTO_TEST_CASE ( DOCI_public_methods ) {
    // Create an AOBasis
    GQCG::Molecule water ("../tests/data/h2o.xyz");
    auto ao_basis = std::make_shared<GQCG::AOBasis>(water, "STO-3G");


    // Create random HamiltonianParameters from One- and TwoElectronOperators (and a transformation matrix) with compatible dimensions
    size_t K = ao_basis->get_number_of_basis_functions();
    GQCG::OneElectronOperator S (Eigen::MatrixXd::Random(K, K));
    GQCG::OneElectronOperator H_core (Eigen::MatrixXd::Random(K, K));
    Eigen::Tensor<double, 4> g_tensor (K, K, K, K);
    g_tensor.setRandom();
    GQCG::TwoElectronOperator g (g_tensor);
    Eigen::MatrixXd C = Eigen::MatrixXd::Random(K, K);
    GQCG::HamiltonianParameters random_hamiltonian_parameters (ao_basis, S, H_core, g, C);

    // Create a compatible Fock space
    GQCG::FockSpace fock_space (K, 3);

    // Create DOCI module
    GQCG::DOCI random_doci (fock_space);

    // Test the public DOCI methods
    Eigen::VectorXd x = random_doci.calculateDiagonal(random_hamiltonian_parameters);
    BOOST_CHECK_NO_THROW(random_doci.constructHamiltonian(random_hamiltonian_parameters));
    BOOST_CHECK_NO_THROW(random_doci.matrixVectorProduct(random_hamiltonian_parameters, x, x));

    // Create an incompatible Fock space
    GQCG::FockSpace fock_space_i (K+1, 3);

    // Create DOCI module
    GQCG::DOCI random_doci_i (fock_space_i);
    BOOST_CHECK_THROW(random_doci_i.constructHamiltonian(random_hamiltonian_parameters), std::invalid_argument);
    BOOST_CHECK_THROW(random_doci_i.matrixVectorProduct(random_hamiltonian_parameters, x, x), std::invalid_argument);
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
   GQCG::DOCI doci (fock_space);

   // Retrieve the Hamiltonian matrix
   Eigen::MatrixXd ham = doci.constructHamiltonian(ham_par);

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
    GQCG::DOCI doci (fock_space);

    // Retrieve the Hamiltonian matrix
    Eigen::MatrixXd ham = doci.constructHamiltonian(ham_par);

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
    GQCG::DOCI doci (fock_space);

    // Retrieve the Hamiltonian matrix
    Eigen::MatrixXd ham = doci.constructHamiltonian(ham_par);

    // Retrieve the eigenvalues
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen_solver (ham);
    auto doci_eigenvalue = eigen_solver.eigenvalues()(0);

    // Calculate the total energy
    double internuclear_repulsion_energy =  3.0036546888874875e+00;  // this comes straight out of the FCIDUMP file
    double test_doci_energy = doci_eigenvalue + internuclear_repulsion_energy;

    BOOST_CHECK(std::abs(test_doci_energy - (reference_doci_energy)) < 1.0e-9);
}


BOOST_AUTO_TEST_CASE ( DOCI_h2o_sto3g_klaas_Davidson ) {

    // Klaas' reference DOCI energy for H2O@STO-3G
    double reference_doci_energy = -74.9671366903;


    // Do a DOCI calculation based on a given FCIDUMP file
    // Create the Hamiltonian Parameters
    auto mol_ham_par = GQCG::readFCIDUMPFile("../tests/data/h2o_sto3g_klaas.FCIDUMP");

    // The species contains 10 electrons and 7 basis functions, this requires a single Fock Space of 7 orbitals and 5 electrons
    GQCG::FockSpace fock_space (7, 5);  // dim = 21

    // Create the DOCI module
    GQCG::DOCI doci (fock_space);
    Eigen::VectorXd diagonal = doci.calculateDiagonal(mol_ham_par);

    // Davidson Solver requires us to specify the macvec:
    numopt::VectorFunction matrixVectorProduct = [&doci, &mol_ham_par, &diagonal](const Eigen::VectorXd& x) { return doci.matrixVectorProduct(mol_ham_par, x, diagonal); };

    // Davidson Solver requires initial guess
    Eigen::VectorXd initial_guess = Eigen::VectorXd::Zero(21);
    initial_guess(0) = 1;

    // Davidson Solver requires diagonal of the DOCI Hamiltonian
    numopt::eigenproblem::DavidsonSolver davidson_solver (matrixVectorProduct, diagonal, initial_guess);
    davidson_solver.solve();

    auto doci_eigenvalue = davidson_solver.get_lowest_eigenvalue();

    // Calculate the total energy
    double internuclear_repulsion_energy =  9.7794061444134091E+00;  // this comes straight out of the FCIDUMP file
    double test_doci_energy = doci_eigenvalue + internuclear_repulsion_energy;

    BOOST_CHECK(std::abs(test_doci_energy - (reference_doci_energy)) < 1.0e-9);
}
