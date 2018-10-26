/**
 *  An executable that calculates the OO-DOCI energy, starting from the RHF orbitals
 */
#include <fstream>
#include <iomanip>
#include <iostream>

#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>

namespace po = boost::program_options;

#include "HamiltonianParameters/HamiltonianParameters_constructors.hpp"
#include "RHF/DIISRHFSCFSolver.hpp"
#include "CISolver/CISolver.hpp"
#include "DOCINewtonOrbitalOptimizer.hpp"



int main (int argc, char** argv) {

    // Input processing
    std::string input_xyz_file;
    std::string basisset;

    po::variables_map variables_map;
    try {
        po::options_description desc ("Options");
        desc.add_options()
        ("help,h", "print help messages")
        ("input,f", po::value<std::string>(&input_xyz_file)->required(), "filename of the .xyz-file")
        ("basis,b", po::value<std::string>(&basisset), "name of the basis set");

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
    // Prepare molecular Hamiltonian parameters in the RHF basis
    GQCP::Molecule molecule (input_xyz_file);
    auto ao_basis = std::make_shared<GQCP::AOBasis>(molecule, basisset);
    auto ao_mol_ham_par = GQCP::constructMolecularHamiltonianParameters(ao_basis);

    GQCP::DIISRHFSCFSolver diis_scf_solver (ao_mol_ham_par, molecule);
    diis_scf_solver.solve();
    auto rhf = diis_scf_solver.get_solution();

    auto mol_ham_par = GQCP::HamiltonianParameters(ao_mol_ham_par, rhf.get_C());


    // Do the DOCI orbital optimization using specified solver options
    auto K = ao_basis->get_number_of_basis_functions();
    GQCP::FockSpace fock_space (K, molecule.get_N()/2);  // provide the number of electron pairs
    GQCP::DOCI doci (fock_space);

    GQCP::DOCINewtonOrbitalOptimizer orbital_optimizer (doci, mol_ham_par);
    numopt::eigenproblem::DenseSolverOptions solver_options;
    orbital_optimizer.solve(solver_options);

    double OO_DOCI_electronic_energy = orbital_optimizer.get_eigenpair().get_eigenvalue();
    double internuclear_repulsion_energy = molecule.calculateInternuclearRepulsionEnergy();


    // Print the energy to an output file
    // Create and open a file: filename.xyz -> filename_doci_rhf_basisset.output
    std::string output_filename = input_xyz_file;
    boost::replace_last(output_filename, ".xyz", std::string("_doci_rhf_") + basisset + std::string(".output"));

    std::ofstream output_file;
    output_file.open(output_filename, std::fstream::out);

    output_file << std::setprecision(15) << OO_DOCI_electronic_energy + internuclear_repulsion_energy << std::endl;

    output_file.close();
}
