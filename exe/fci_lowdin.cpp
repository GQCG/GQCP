/**
 *  An executable that calculates the FCI energy
 */
#include <fstream>
#include <iomanip>
#include <iostream>

#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>

namespace po = boost::program_options;

#include "HamiltonianParameters/HamiltonianParameters_constructors.hpp"
#include "CISolver/CISolver.hpp"
#include "HamiltonianBuilder/FCI.hpp"
#include "FockSpace/ProductFockSpace.hpp"
#include "RDM/RDMCalculator.hpp"



int main (int argc, char** argv) {

    // Input processing
    std::string input_xyz_file;
    std::string basisset;
    size_t N_alpha;
    size_t N_beta;

    po::variables_map variables_map;
    try {
        po::options_description desc ("Options");
        desc.add_options()
        ("help,h", "print help messages")
        ("input,f", po::value<std::string>(&input_xyz_file)->required(), "filename of the .xyz-file")
        ("N_alpha,a", po::value<size_t>(&N_alpha)->required(), "number of alpha electrons")
        ("N_beta,b", po::value<size_t>(&N_beta)->required(), "number of beta electrons")
        ("basis,s", po::value<std::string>(&basisset)->required(), "name of the basis set");


        po::store(po::parse_command_line(argc, argv, desc), variables_map);

        if (variables_map.count("help")) {
            std::cout << "ONV expansion 1- and 2-RDMs" << std::endl << desc << std::endl;
            std::exit(0);
        }

        po::notify(variables_map);
    } catch (po::error& e) {
        std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
        return 1;
    } catch(...) {
        std::cerr << "ERROR: you have not specified all arguments. Please use the -h flag for more information." << std::endl << std::endl;
    }


    // Actual calculations
    // Prepare molecular Hamiltonian parameters in the Löwdin basis
    GQCP::Molecule molecule (input_xyz_file);
    auto ao_basis = std::make_shared<GQCP::AOBasis>(molecule, basisset);
    auto mol_ham_par = GQCP::constructMolecularHamiltonianParameters(ao_basis);  // in the AO basis
    mol_ham_par.LowdinOrthonormalize();  // now in the Löwdin basis


    // Solve the FCI eigenvalue problem using the dense algorithm
    auto K = ao_basis->get_number_of_basis_functions();
    GQCP::ProductFockSpace fock_space (K, N_alpha, N_beta);
    GQCP::FCI fci (fock_space);
    GQCP::CISolver ci_solver (fci, mol_ham_par);
    numopt::eigenproblem::DenseSolverOptions dense_solver_options;
    ci_solver.solve(dense_solver_options);

    auto fci_energy = ci_solver.get_eigenpair().get_eigenvalue();
    double internuclear_repulsion_energy = molecule.calculateInternuclearRepulsionEnergy();


    // Calculate the 1-RDM in the NO basis
    auto fci_coefficients = ci_solver.get_eigenpair().get_eigenvector();
    GQCP::RDMCalculator fci_rdm (fock_space);
    GQCP::OneRDM D = fci_rdm.calculate1RDMs(fci_coefficients).one_rdm;
    D.diagonalize();


    // Print the energy to an output file
    // Create and open a file: filename.xyz -> filename_doci_rhf_basisset.output
    std::string output_filename = input_xyz_file;
    boost::replace_last(output_filename, ".xyz", std::string("_fci_lowdin_") + basisset + std::string(".output"));

    std::ofstream output_file;
    output_file.open(output_filename, std::fstream::out);

    output_file << "TOTAL ENERGY: " << std::setprecision(15) << fci_energy + internuclear_repulsion_energy << std::endl << std::endl;

    output_file << "1-DM IN NATURAL ORBITAL BASIS\n" << D.get_matrix_representation().diagonal() << std::endl << std::endl;

    output_file.close();
}
