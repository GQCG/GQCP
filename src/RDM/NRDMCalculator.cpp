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
//    ONV ket = fock_space.get_ONV(0);



    size_t I = 0;
    while (I < dim) {  // loop over all bra addresses

        int sign = 1;  // (re)set the total sign


        // Annihilate the bra on the bra indices
        size_t i = 0;
        while (i < bra_indices.size()) {

            size_t bra_index = bra_indices[i];
            if (!bra.annihilate(bra_index, sign)) {

                // If we can't annihilate on the current bra, we should immediately go to the next bra
                if (I < dim-1) {  // prevent the last permutation to occur
                    fock_space.setNext(bra);
                    I++;
                }
                i = 0;  // for every new I, we have to restart the bra indices as well
            } else {
                std::cout << "We annihilated on " << bra_index << std::endl;
                std::cout << "The current bra is now: " << bra.asString() << std::endl;
            }
        }

        std::cout << bra.asString() << std::endl;
    }






//    for (size_t I = 0; I < dim; I++) {  // loop over all bra addresses
//
//
//
//
//
////        for (size_t J = 0; J < dim; J++) {  // loop over all ket addresses
////
////            // Annihilate the ket on the ket indices
////            for (const auto& ket_index : ket_indices) {
////
////                if (!ket.annihilate(ket_index, sign)) {
////
////                    // If we can't annihilate on the current ket, we should immediately go to the next bra
////                    if (J < dim-1) {  // prevent the last permutation to occur
////                        fock_space.setNext(ket);
////                        J++;  // we aren't at the end of the for-loop, so increment J manually
////                    }
////                    continue;
////                }
////            }
////
////            // Right here, we're sure all the annihilations on the bra and the ket have been successful
////            if (bra.countNumberOfDifferences(ket) == 0) {
////                value += sign * coeff(I) * coeff(J);
////            }
////
////
////
////
////            if (J < dim-1) {  // prevent the last permutation to occur
////                fock_space.setNext(ket);
////            }
////
////        }  // loop over J
//
//
//
//        if (I < dim-1) {  // prevent the last permutation to occur
//            fock_space.setNext(bra);
//        }
//    }  // loop over I


    return value;
}





}  // namespace GQCP
