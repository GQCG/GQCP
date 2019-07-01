
#include "HamiltonianParameters/HamiltonianParameters.hpp"
#include "CISolver/CISolver.hpp"
#include "HamiltonianBuilder/FCI.hpp"
#include "FockSpace/ProductFockSpace.hpp"
#include "RDM/RDMCalculator.hpp"

namespace GQCP {
    class FullConfigurationInteractionDriver {
    private:
        std::shared_ptr<GQCP::Molecule> molecule;
        std::shared_ptr<GQCP::CISolver> solver;

    public:
        FullConfigurationInteractionDriver(std::string xyz_filename, std::string basis_set, size_t num_alpha, size_t num_beta);

        double get_energy();
    };
}
