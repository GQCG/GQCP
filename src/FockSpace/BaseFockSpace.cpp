#include "FockSpace/BaseFockSpace.hpp"


namespace GQCG {

/*
 * PROTECTED CONSTRUCTORS
 */

/*
 *  Protected constructor given a @param K
 */
BaseFockSpace::BaseFockSpace(size_t K) :
    K(K)
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
 *  Creates a Hartree-Fock guess
 */
Eigen::VectorXd BaseFockSpace::hartreeFockGuess(){
    Eigen::VectorXd guess = Eigen::VectorXd::Zero(this->dim);
    guess(0) = 1;
    return guess;
}



}  // namespace GQCG
