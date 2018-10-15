#include "WaveFunction.hpp"


namespace GQCG {


/*
 * CONSTRUCTORS
 */

WaveFunction::WaveFunction(BaseFockSpace& base_fock_space, const Eigen::VectorXd& coefficients) :
    fock_space (&base_fock_space),
    coefficients (coefficients)
{}


}  // namespace GQCG
