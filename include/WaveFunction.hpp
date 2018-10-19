#ifndef GQCP_WAVEFUNCTION_HPP
#define GQCP_WAVEFUNCTION_HPP


#include "FockSpace/BaseFockSpace.hpp"

#include "common.hpp"



namespace GQCP {


/**
 *  WaveFunction contains the expansion coefficients in its given FockSpace
 */
class WaveFunction {
private:
    BaseFockSpace* fock_space;
    Eigen::VectorXd coefficients;  // Expansion coefficients of a wave function in the Fock space

public:
    // CONSTRUCTORS
    WaveFunction(BaseFockSpace& base_fock_space, const Eigen::VectorXd& coefficients);


    // GETTERS
    Eigen::VectorXd get_coefficients() const { return coefficients; }
    BaseFockSpace& get_fock_space() const { return *fock_space; }
};


}  // namespace GQCP


#endif  // GQCP_WAVEFUNCTION_HPP
