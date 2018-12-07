#ifndef NRDMCalculator_hpp
#define NRDMCalculator_hpp


#include "FockSpace/FockSpace.hpp"
#include <iostream>


namespace GQCP {


/**
 *  A class that can be used to calculate elements of N-th order density matrices for full CI wave functions (i.e. a FockSpace with spin orbitals)
 */
class NRDMCalculator {
private:
    FockSpace fock_space;  // the FCI Fock space with spin orbitals


public:
    // CONSTRUCTORS
    /**
     *  @param fock_space       the FCI Fock space with spin orbitals
     */
    explicit NRDMCalculator(const FockSpace& fock_space);


    // PUBLIC METHODS
    /**
     *  @param bra_indices      the indices of the orbitals that should be annihilated on the left (on the bra)
     *  @param ket_indices      the indices of the orbitals that should be annihilated on the right (on the ket)
     *  @param coeff            the expansion coefficient vector
     *
     *  @return an element of the N-RDM, as specified by the given bra and ket indices
     *
     *      calculateElement({0, 1}, {2, 1}) would calculate d^{(2)} (0, 1, 1, 2): the operator string would be a^\dagger_0 a^\dagger_1 a_2 a_1
     */
    double calculateElement(const std::vector<size_t>& bra_indices, const std::vector<size_t>& ket_indices, const Eigen::VectorXd& coeff) const;
};




}  // namespace GQCP

#endif /* NRDMCalculator_hpp */
