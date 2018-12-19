/**
 *  An executable that calculates the Hubbard energy/energies, starting from a triagonal input specifying the hopping matrix. The energy/energies are printed to the console.
 *
 *  Example execution:
 *      ./hubbard -m "comma-separated upper triagonal" -a "number of alpha electrons" -b "number of beta electrons" -K "number of sites"
 *
 *  For example:
 *      ./hubbard -m "1,2,3,4,5,6" -a 1 -b 2 -K 3
 *  Default amount of eigenvalues is 1, this can be changed with the x flag:
 *      ./hubbard -m "1,2,3,4,5,6" -x 2 -a 1 -b 2 -K 3
 */




#include <fstream>
#include <iomanip>
#include <iostream>

#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>

namespace po = boost::program_options;

#include "HamiltonianParameters/HamiltonianParameters.hpp"
#include "CISolver/CISolver.hpp"
#include "HamiltonianBuilder/Hubbard.hpp"


int main (int argc, char** argv) {

    // Input processing
    std::string csline;
    size_t N_alpha;
    size_t N_beta;
    size_t K;
    size_t N_eigenvalues;
    po::variables_map variables_map;

    try {
        po::options_description desc ("Options");
        desc.add_options()
        ("help,h", "print help messages")
        ("csline,m", po::value<std::string>(&csline)->required(), "Comma-separated upper triagonal line")
        ("N_lowest_states,x", po::value<size_t>(&N_eigenvalues)->default_value(1), "Number of lowest states")
        ("N_alpha,a", po::value<size_t>(&N_alpha)->required(), "number of alpha electrons")
        ("N_beta,b", po::value<size_t>(&N_beta)->required(), "number of beta electrons")
        ("K,K", po::value<size_t>(&K)->required(), "number of sites");

        po::store(po::parse_command_line(argc, argv, desc), variables_map);

        if (variables_map.count("help")) {
            std::cout << "Hubbard eigenvalues" << std::endl << desc << std::endl;
            std::exit(0);
        }

        po::notify(variables_map);
    } catch (po::error& e) {
        std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
        return 1;
    } catch(...) {
        std::cerr << "ERROR: you have not specified all arguments. Please use the -h flag for more information." << std::endl << std::endl;
    }

    std::string triagonal_line = csline;

    // Read the upper triagonal of the hopping matrix
    Eigen::VectorXd triagonal;
    if (csline.empty()){
        throw std::invalid_argument("Comma-separated was empty!");
    }

    // Split comma-separated line
    std::vector<std::string> splitted_line;
    boost::split(splitted_line, triagonal_line, boost::is_any_of(","));

    std::vector<double> triagonal_data;
    for (const std::string& x : splitted_line) {
        triagonal_data.push_back(std::stod(x));
    }

    triagonal = Eigen::Map<Eigen::VectorXd>(triagonal_data.data(), triagonal_data.size());
    GQCP::HoppingMatrix H = GQCP::HoppingMatrix::FromUpperTriangle(triagonal);

    // Actual calculations
    auto ham_par = GQCP::HamiltonianParameters::Hubbard(H);
    if (ham_par.get_K() != K) {
        throw std::invalid_argument("The given number of sites does not match the triagonal");
    }


    // Initialize and solve the Hubbard eigenvalue problem
    GQCP::ProductFockSpace fock_space (K, N_alpha, N_beta);
    GQCP::Hubbard hubbard (fock_space);
    GQCP::CISolver solver (hubbard, ham_par);
    GQCP::DenseSolverOptions dense_solver_options;
    dense_solver_options.number_of_requested_eigenpairs = N_eigenvalues;
    solver.solve(dense_solver_options);


    // Print the energy to the console
    std::cout << std::setprecision(15);
    for (const GQCP::Eigenpair& eigenpair : solver.get_eigenpairs()) {
        std::cout << eigenpair.get_eigenvalue() << std::endl;
    }
}
