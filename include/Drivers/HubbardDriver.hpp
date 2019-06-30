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



#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>

#include "HamiltonianParameters/HamiltonianParameters.hpp"
#include "CISolver/CISolver.hpp"
#include "HamiltonianBuilder/Hubbard.hpp"
#include "RDM/RDMCalculator.hpp"

namespace GQCP {
    class HubbardDriver {
    private:
        std::shared_ptr<GQCP::ProductFockSpace> fock_space;
        std::shared_ptr<GQCP::Hubbard> hubbard;
        std::shared_ptr<GQCP::CISolver> solver;
        GQCP::DenseSolverOptions dense_solver_options;

    public:
        HubbardDriver(std::string csline, size_t num_states, size_t num_alpha, size_t num_beta, size_t K);

        std::vector<double> get_energies();

        std::vector<Eigen::MatrixXd> get_first_order_rdms();
    };
}
