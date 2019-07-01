#include "Drivers/HubbardDriver.hpp"

namespace GQCP {

    HubbardDriver::HubbardDriver(std::string csline, size_t num_states, size_t num_alpha, size_t num_beta, size_t K) {

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
        this->fock_space = std::make_shared<GQCP::ProductFockSpace>(K, num_alpha, num_beta);
        this->hubbard = std::make_shared<GQCP::Hubbard>(*(this->fock_space));
        this->solver = std::make_shared<GQCP::CISolver>(*(this->hubbard), ham_par);
        this->dense_solver_options.number_of_requested_eigenpairs = num_states;
        this->solver->solve(this->dense_solver_options);
    }

    std::vector<double> GQCP::HubbardDriver::get_energies() {
        std::vector<double> energies;
        for (const GQCP::Eigenpair &eigenpair : this->solver->get_eigenpairs()) {
            energies.push_back(eigenpair.get_eigenvalue());
        }
        return energies;
    }

    std::vector<Eigen::MatrixXd> GQCP::HubbardDriver::get_first_order_rdms() {
        std::vector<Eigen::MatrixXd> rdms;
        for (const GQCP::Eigenpair &eigenpair : this->solver->get_eigenpairs()) {
            GQCP::RDMCalculator rdm_calculator(*(this->fock_space));
            rdm_calculator.set_coefficients(eigenpair.get_eigenvector());
            rdms.push_back(static_cast<Eigen::MatrixXd>(rdm_calculator.calculate1RDMs().one_rdm));
        }
        return rdms;
    }
}
