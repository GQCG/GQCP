#ifndef GQCG_WAVEFUNCTION_HPP
#define GQCG_WAVEFUNCTION_HPP


#include "FockSpace/BaseFockSpace.hpp"

#include "common.hpp"



namespace GQCG {


class WaveFunction {
private:
    BaseFockSpace* base_fock_space;
    Eigen::VectorXd coefficients;  // Expansion coefficients of a wave function in the Fock space

public:
    // CONSTRUCTORS
    WaveFunction(BaseFockSpace& base_fock_space, const Eigen::VectorXd& coefficients) : base_fock_space(&base_fock_space), coefficients(coefficients){}


    // GETTERS
    Eigen::VectorXd get_coefficients() { return coefficients; }
    BaseFockSpace& get_fock_space() { return *base_fock_space; }
};


}  // namespace GQCG
#endif  // GQCG_WAVEFUNCTION_HPP
