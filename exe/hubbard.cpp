/**
 *  An executable that calculates the Hubbard energy, starting from a triagonal input specifying the hopping matrix
 *
 *  Example execution:
 *      ./hubbard -f "filename" -a "number of alpha electrons" -b "number of beta electrons" -K "number of sites"
 *
 *  For example:
 *      ./hubbard -f hubbard_in -a 1 -b 2 -K 3
 *  where hubbard_in is a comma-separated file:
 *      for K=3: 1,2,3,4,5,6
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
    size_t N_alpha;
    size_t N_beta;
    size_t K;

    po::variables_map variables_map;

    try {
        po::options_description desc ("Options");
        desc.add_options()
        ("help,h", "print help messages")
        ("input,f", po::value<std::string>(&hubbard_input_file)->required(), "filename of the input-file")
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


    // Actual calculations
    // Read the upper triagonal of the hopping matrix

    std::ifstream file (hubbard_input_file);
    Eigen::VectorXd triagonal;
    if (file.is_open()) {

        // Split comma-separated line
        std::string triagonal_line;
        std::getline(file, triagonal_line);
        std::vector<std::string> splitted_line;
        boost::split(splitted_line, triagonal_line, boost::is_any_of(","));

        std::vector<double> triagonal_data;
        for (const std::string& x : splitted_line) {
            triagonal_data.push_back(std::stod(x));
        }

        triagonal = Eigen::Map<Eigen::VectorXd>(triagonal_data.data(), triagonal_data.size());

        file.close();
    } else {
        throw std::invalid_argument("File was not found!");
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
    solver.solve(dense_solver_options);


    // Print the energy to an output file: hubbard_input_file.out
    std::string output_filename = hubbard_input_file + ".out";
    std::ofstream output_file;
    output_file.open(output_filename, std::fstream::out);

    output_file << std::setprecision(15) << solver.get_eigenpair().get_eigenvalue() << std::endl;

    output_file.close();
}
