#include "RDM/NRDMCalculator.hpp"



#include <iostream>



namespace GQCP {


/*
 *  CONSTRUCTORS
 */
/**
 *  @param fock_space       the FCI Fock space with spin orbitals
 */
NRDMCalculator::NRDMCalculator(const FockSpace& fock_space) :
    fock_space (fock_space)
{}




/*
 *  PUBLIC METHODS
 */
/**
 *  @param bra_indices      the indices of the orbitals that should be annihilated on the left (on the bra)
 *  @param ket_indices      the indices of the orbitals that should be annihilated on the right (on the ket)
 *  @param coeff            the expansion coefficient vector
 */
double NRDMCalculator::calculateElement(const std::vector<size_t>& bra_indices, const std::vector<size_t>& ket_indices, const Eigen::VectorXd& coeff) const {


    std::cout << "bra indices:" << std::endl;
    for (const auto& bra_index : bra_indices) {
        std::cout << bra_index << std::endl;
    }
    std::cout << std::endl;

    std::cout << "ket indices" << std::endl;
    for (const auto& ket_index : ket_indices) {
        std::cout << ket_index << std::endl;
    }



    double value = 0.0;
    FockSpace fock_space = this->fock_space;  // make a copy because this method is marked const
    size_t dim = fock_space.get_dimension();


    ONV bra = fock_space.get_ONV(0);
    std::cout << "Initial bra: " << bra.asString() << std::endl;
    ONV ket = fock_space.get_ONV(0);
    std::cout << "Initial ket: " << ket.asString() << std::endl;



    size_t I = 0;
    while (I < dim) {  // loop over all bra addresses

        int sign = 1;  // (re)set the total sign


        // Annihilate the bra on the bra indices
        if (!bra.annihilateAll(bra_indices, sign)) {

            // Go to the beginning of the outer while loop with the next bra
            if (I < dim-1) {  // prevent the last permutation to occur
                fock_space.setNext(bra);
                I++;
                continue;
            }
        }

        std::cout << "The fully-annihilated bra is now: " << bra.asString() << std::endl;


        size_t J = 0;
        while (J < dim) {

            // Annihilate the ket on the ket indices
            if (!ket.annihilateAll(ket_indices, sign)) {

                // Go to the beginning of this (the inner) while loop with the next bra
                if (J < dim-1) {  // prevent the last permutation to occur
                    fock_space.setNext(ket);
                    J++;
                    continue;
                }
            }


            std::cout << "The fully-annihilated ket is now: " << ket.asString() << std::endl;

            if (bra == ket) {
                value += sign * coeff(I) * coeff(J);
            }

        }  // while J loop

    }  // while I loop


    return value;
}





}  // namespace GQCP
