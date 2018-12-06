#ifndef NRDMCalculator_hpp
#define NRDMCalculator_hpp


#include "FockSpace/FockSpace.hpp"



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
     */
    double calculateElement(const std::vector<size_t>& bra_indices, const std::vector<size_t>& ket_indices, const Eigen::VectorXd& coeff) const;


    // provide an operator() to simplify the API
};




}  // namespace GQCP

#endif /* NRDMCalculator_hpp */
