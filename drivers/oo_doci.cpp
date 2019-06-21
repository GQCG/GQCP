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


#include "HamiltonianParameters/HamiltonianParameters.hpp"
#include "RHF/DIISRHFSCFSolver.hpp"
#include "CISolver/CISolver.hpp"
#include "OrbitalOptimization/DOCINewtonOrbitalOptimizer.hpp"
#include "OrbitalOptimization/Localization/ERNewtonLocalizer.hpp"
#include "OrbitalOptimization/Localization/ERJacobiLocalizer.hpp"
#include "math/optimization/IterativeIdentitiesHessianModifier.hpp"


/**
 *  The main function for the OO-DOCI executable
 */
int main (int argc, char** argv) {

    // Input processing
    std::string input_xyz_file;
    std::string basisset;
    bool user_wants_localization = false;
    bool user_wants_davidson = false;

    namespace po = boost::program_options;
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
    auto molecule = GQCP::Molecule::Readxyz(input_xyz_file);
    size_t N_P = molecule.get_N()/2;
    output_file << "Molecule geometry" << std::endl;
    output_file << molecule << std::endl;


    output_file << "Basisset: " << basisset << std::endl << std::endl;
    auto ao_mol_ham_par = GQCP::HamiltonianParameters<double>::Molecular(molecule, basisset);
    size_t K = ao_mol_ham_par.get_K();


    GQCP::DIISRHFSCFSolver diis_scf_solver (ao_mol_ham_par, molecule);
    diis_scf_solver.solve();
    auto rhf = diis_scf_solver.get_solution();

    output_file << "Transformation matrix to the RHF orbitals: " << std::endl << rhf.get_C() << std::endl << std::endl;
    auto mol_ham_par = GQCP::HamiltonianParameters<double>(ao_mol_ham_par, rhf.get_C());


    // Localize the RHF-orbitals
    auto hessian_modifier = std::make_shared<GQCP::IterativeIdentitiesHessianModifier>();
    if (user_wants_localization) {

        // Newton to get to the first local minimum
        GQCP::ERNewtonLocalizer first_newton_localizer (N_P, hessian_modifier);
        first_newton_localizer.optimize(mol_ham_par);


        // Check if Jacobi finds another minimum
        GQCP::ERJacobiLocalizer jacobi_localizer (N_P);
        auto optimal_jacobi_with_scalar = jacobi_localizer.calculateOptimalJacobiParameters(mol_ham_par);
        if (optimal_jacobi_with_scalar.second > 0) {  // if a Jacobi rotation can find an increase, do it
            const auto U = GQCP::SquareMatrix<double>::FromJacobi(optimal_jacobi_with_scalar.first, mol_ham_par.get_K());
            mol_ham_par.rotate(U);
        }


        // Newton to get to the next local minimum
        GQCP::ERNewtonLocalizer second_newton_localizer (N_P, hessian_modifier);
        first_newton_localizer.optimize(mol_ham_par);


        output_file << "Total transformation matrix to the ER-localized orbitals: " << std::endl << mol_ham_par.get_T_total() << std::endl << std::endl;
    }


    // Do the DOCI orbital optimization
    GQCP::FockSpace fock_space (K, N_P);
    GQCP::DOCI doci (fock_space);

    std::shared_ptr<GQCP::BaseSolverOptions> solver_options;
    if (user_wants_davidson) {
        solver_options = std::make_shared<GQCP::DavidsonSolverOptions>(fock_space.HartreeFockExpansion());
    } else {
        solver_options = std::make_shared<GQCP::DenseSolverOptions>();
    }

    GQCP::DOCINewtonOrbitalOptimizer orbital_optimizer (doci, *solver_options, hessian_modifier);
    orbital_optimizer.optimize(mol_ham_par);
    double OO_DOCI_electronic_energy = orbital_optimizer.get_eigenpair().get_eigenvalue();


    output_file << "Total transformation matrix to the OO-DOCI orbitals: " << std::endl << mol_ham_par.get_T_total() << std::endl << std::endl;

    output_file << "Total OO-DOCI energy (internuclear repulsion energy added): " << std::setprecision(15) << OO_DOCI_electronic_energy + mol_ham_par.get_scalar() << std::endl;

    output_file.close();
}
