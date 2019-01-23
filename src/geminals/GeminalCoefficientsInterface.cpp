#include "geminals/GeminalCoefficientsInterface.hpp"


namespace GQCP {


/*
 *  DESTRUCTOR
 */
GeminalCoefficientsInterface::~GeminalCoefficientsInterface() {}


/*
 *  PUBLIC METHODS
 */

/**
 *  @param fock_space       the seniority-zero Fock space the wave function should live in
 *
 *  @return the wave function expansion corresponding to the geminal coefficients
 */
WaveFunction GeminalCoefficientsInterface::toWaveFunction(const FockSpace& fock_space) const {

    // The FockSpace can't be marked const, as makeONV() and setNext() are non-const methods

    Eigen::VectorXd coefficients = Eigen::VectorXd::Zero(fock_space.get_dimension());  // coefficient vector
    ONV onv = fock_space.makeONV(0);  // start with address 0
    for (size_t I = 0; I < fock_space.get_dimension(); I++) {

        coefficients(I) = this->overlap(onv);

        if (I < fock_space.get_dimension() - 1) {  // skip the last permutation
            fock_space.setNextONV(onv);
        }
    }

    return WaveFunction(fock_space, coefficients);
}


}  // namespace GQCP
