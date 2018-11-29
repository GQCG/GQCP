/**
 *  An executable that calculates the OO-DOCI energy, starting from the RHF orbitals
 *  Optionally:
 *      - the RHF orbitals can be localized using Edmiston-Ruedenberg localization before starting the OO-DOCI procedure
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
#include "Localization/ERNewtonLocalizer.hpp"
#include "Localization/ERJacobiLocalizer.hpp"


/**
 *  The main function for the OO-DOCI executable
 */
int main (int argc, char** argv) {

    // Input processing
    std::string input_xyz_file;
    std::string basisset;
    bool user_wants_localization = false;
    bool user_wants_davidson = false;

    po::variables_map variables_map;
    try {
        po::options_description desc ("Options");
        desc.add_options()
        ("help,h", "print help messages")
        ("input,f", po::value<std::string>(&input_xyz_file)->required(), "filename of the .xyz-file")
        ("basis,b", po::value<std::string>(&basisset)->required(), "name of the basis set")
        ("solver,s", po::bool_switch(&user_wants_davidson)->required(), "if the Davidson diagonalization algorithm should be used")
        ("localize,l", po::bool_switch(&user_wants_localization), "if the RHF orbitals should be localized using the ER-localization");

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


    // Create and open a file: filename.xyz -> filename.output
    std::string output_filename = input_xyz_file;
    boost::replace_last(output_filename, ".xyz", ".output");

    std::ofstream output_file;
    output_file.open(output_filename, std::fstream::out);



    // Actual calculations
    // Prepare molecular Hamiltonian parameters in the RHF basis
    GQCP::Molecule molecule (input_xyz_file);
    double internuclear_repulsion_energy = molecule.calculateInternuclearRepulsionEnergy();
    size_t N_P = molecule.get_N()/2;
    output_file << "Molecule geometry" << std::endl;
    output_file << molecule << std::endl;


    auto ao_basis = std::make_shared<GQCP::AOBasis>(molecule, basisset);
    output_file << "Basisset: " << basisset << std::endl << std::endl;
    size_t K = ao_basis->get_number_of_basis_functions();

    auto ao_mol_ham_par = GQCP::constructMolecularHamiltonianParameters(ao_basis);


    GQCP::DIISRHFSCFSolver diis_scf_solver (ao_mol_ham_par, molecule);
    diis_scf_solver.solve();
    auto rhf = diis_scf_solver.get_solution();

    output_file << "Transformation matrix to the RHF orbitals: " << std::endl << rhf.get_C() << std::endl << std::endl;
    auto mol_ham_par = GQCP::HamiltonianParameters(ao_mol_ham_par, rhf.get_C());


    // Localize the RHF-orbitals
    if (user_wants_localization) {

        // Newton to get to the first local minimum
        GQCP::ERNewtonLocalizer first_newton_localizer (N_P);
        first_newton_localizer.localize(mol_ham_par);


        // Check if Jacobi finds another minimum (1 iteration)
        GQCP::ERJacobiLocalizer jacobi_localizer (N_P, 1.0, 1);  // Do 1 iteration, regardless of threshold
        jacobi_localizer.localize(mol_ham_par);


        // Newton to get to the next local minimum
        GQCP::ERNewtonLocalizer second_newton_localizer (N_P);
        first_newton_localizer.localize(mol_ham_par);


        output_file << "Total transformation matrix to the ER-localized orbitals: " << std::endl << mol_ham_par.get_C() << std::endl << std::endl;
    }


    // Do the DOCI orbital optimization
    GQCP::FockSpace fock_space (K, N_P);
    GQCP::DOCI doci (fock_space);

    GQCP::DOCINewtonOrbitalOptimizer orbital_optimizer (doci, mol_ham_par);

    if (user_wants_davidson) {
        numopt::eigenproblem::DavidsonSolverOptions solver_options (fock_space.HartreeFockExpansion());
        orbital_optimizer.solve(solver_options);
    } else {
        numopt::eigenproblem::DenseSolverOptions solver_options;
        orbital_optimizer.solve(solver_options);
    }
    output_file << "Total transformation matrix to the OO-DOCI orbitals: " << std::endl << mol_ham_par.get_C() << std::endl << std::endl;

    double OO_DOCI_electronic_energy = orbital_optimizer.get_eigenpair().get_eigenvalue();
    output_file << "Total OO-DOCI energy (internuclear repulsion energy added): " << std::setprecision(15) << OO_DOCI_electronic_energy + internuclear_repulsion_energy << std::endl;

    output_file.close();
}
