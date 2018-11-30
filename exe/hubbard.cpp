/**
 *  An executable that calculates the Hubbard energy(-ies), starting from a triagonal input specifying the hopping matrix
 *
 *  Example execution:
 *      ./hubbard -f "filename" -a "number of alpha electrons" -b "number of beta electrons" -K "number of sites"
 *
 *  For example:
 *      ./hubbard -f hubbard_in -a 1 -b 2 -K 3
 *  where hubbard_in is a comma-separated file:
 *      for K=3: 1,2,3,4,5,6
 *  Default amount of eigenvalues is 1, this can be changed with the x flag:
 *      ./hubbard -f hubbard_in -x 2 -a 1 -b 2 -K 3
 *
 *
 *  Alternatively one can give a non-existing filename and add a comma-separated line to the program arguments:
 *      ./hubbard -f hubbard_out -x 2 -m "1,2,3,4,5,6" -a 1 -b 2 -K 3
 */



#include <fstream>
#include <iomanip>
#include <iostream>

#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>

namespace po = boost::program_options;

#include "HamiltonianParameters/HamiltonianParameters_constructors.hpp"
#include "CISolver/CISolver.hpp"
#include "HamiltonianBuilder/Hubbard.hpp"
#include <string>


int main (int argc, char** argv) {

    // Input processing
    std::string hubbard_input_file;
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
        ("input,f", po::value<std::string>(&hubbard_input_file)->required(), "filename of the output-file and input-file (if no csline is given)")
        ("csline,m", po::value<std::string>(&csline)->default_value(""), "Comma separated upper triagonal line replacing input-file")
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

    bool read_from_file = false;
    std::string triagonal_line;
    std::ifstream file;
    // Actual calculations
    // Read the upper triagonal of the hopping matrix
    Eigen::VectorXd triagonal;
    if (!csline.empty()){
        triagonal_line = csline;
    } else {
        file.open(hubbard_input_file);
        if (file.is_open()) {
            read_from_file = true;
            std::getline(file, triagonal_line);
        } else {
            throw std::invalid_argument("No input-file was found, while no csline was given.");
        }
    }

    // Split comma-separated line
    std::vector<std::string> splitted_line;
    boost::split(splitted_line, triagonal_line, boost::is_any_of(","));

    std::vector<double> triagonal_data;
    for (const std::string& x : splitted_line) {
        triagonal_data.push_back(std::stod(x));
    }

    triagonal = Eigen::Map<Eigen::VectorXd>(triagonal_data.data(), triagonal_data.size());

    if( read_from_file ){
        file.close();
    }

    auto ham_par = GQCP::constructHubbardParameters(triagonal);
    if (ham_par.get_K() != K) {
        throw std::invalid_argument("The given number of sites does not match the triagonal");
    }

    // Initialize and solve the Hubbard eigenvalue problem
    GQCP::ProductFockSpace fock_space (K, N_alpha, N_beta);
    GQCP::Hubbard hubbard (fock_space);
    GQCP::CISolver solver (hubbard, ham_par);
    numopt::eigenproblem::DenseSolverOptions dense_solver_options;
    dense_solver_options.number_of_requested_eigenpairs = N_eigenvalues;
    solver.solve(dense_solver_options);


    // Print the energy to an output file: hubbard_input_file.out
    std::string output_filename = hubbard_input_file + ".out";
    std::ofstream output_file;
    output_file.open(output_filename, std::fstream::out);

    output_file << std::setprecision(15);
    for (const numopt::eigenproblem::Eigenpair& eigenpair : solver.get_eigenpairs()) {
        output_file << eigenpair.get_eigenvalue() << std::endl;
    }

    output_file.close();
}
