#include "Drivers/HubbardDriver.hpp"

namespace GQCP {

    HubbardDriver::HubbardDriver(std::string csline, size_t num_states, size_t num_alpha, size_t num_beta, size_t K) :
    csline (csline),
    num_states(num_states),
    num_alpha(num_alpha),
    num_beta(num_beta),
    K(K) {

        std::string triagonal_line = csline;

        // Read the upper triagonal of the hopping matrix
        GQCP::VectorX<double> triagonal;
        if (csline.empty()) {
            throw std::invalid_argument("hubbard(driver): Comma-separated was empty!");
        }

        // Split comma-separated line
        std::vector<std::string> splitted_line;
        boost::split(splitted_line, triagonal_line, boost::is_any_of(","));

        std::vector<double> triagonal_data;
        for (const std::string &x : splitted_line) {
            triagonal_data.push_back(std::stod(x));
        }

        triagonal = Eigen::Map<Eigen::VectorXd>(triagonal_data.data(), triagonal_data.size());
        GQCP::HoppingMatrix H = GQCP::HoppingMatrix::FromUpperTriangle(triagonal);

        // Actual calculations
        auto ham_par = GQCP::HamiltonianParameters<double>::Hubbard(H);
        if (ham_par.get_K() != K) {
            throw std::invalid_argument("hubbard(driver): The given number of sites does not match the triagonal");
        }


        // Initialize and solve the Hubbard eigenvalue problem
        GQCP::ProductFockSpace fock_space(K, num_alpha, num_beta);
        GQCP::Hubbard hubbard(fock_space);
        GQCP::CISolver solver(hubbard, ham_par);
        GQCP::DenseSolverOptions dense_solver_options;
        dense_solver_options.number_of_requested_eigenpairs = num_states;
        solver.solve(dense_solver_options);
        for (const GQCP::Eigenpair &eigenpair : solver.get_eigenpairs()) {
            this->energies.push_back(eigenpair.get_eigenvalue());
        }
    }
}