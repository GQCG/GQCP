#include "FockSpace/BaseFockSpace.hpp"


namespace GQCG {

/*
 * PROTECTED CONSTRUCTORS
 */

/*
 *  Protected constructor given a @param K
 */
BaseFockSpace::BaseFockSpace(size_t K, size_t dim) :
    K(K),
    dim(dim)
{}



/*
 *  DESTRUCTOR
 */

/**
 *  Provide a pure virtual destructor to make the class abstract
 */
BaseFockSpace::~BaseFockSpace() {}



/*
 *  PUBLIC
 */

/**
 *  Creates a Hartree-Fock coefficient expansion (single Slater expansion of the first configuration in the Fock space)
 */
Eigen::VectorXd BaseFockSpace::HartreeFockExpansion() {
    Eigen::VectorXd expansion = Eigen::VectorXd::Zero(this->dim);
    expansion(0) = 1;  // first configuration is position 0 (conventional ordering of the Fock space)
    return expansion;
}



}  // namespace GQCG
