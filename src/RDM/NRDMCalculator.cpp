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
 *
 *  @return an element of the N-RDM, as specified by the given bra and ket indices
 *
 *      calculateElement({0, 1}, {2, 1}) would calculate d^{(2)} (0, 1, 1, 2)
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

    // The ket indices should be reversed because the annihilators on the ket should be read from right to left
    std::vector<size_t> ket_indices_reversed = ket_indices;
    std::reverse(ket_indices_reversed.begin(), ket_indices_reversed.end());



    double value = 0.0;
    int sign = 1;
    FockSpace fock_space = this->fock_space;  // make a copy because this method is marked const
    size_t dim = fock_space.get_dimension();


    ONV bra = fock_space.get_ONV(0);
    size_t I = 0;
    while (I < dim) {  // loop over all bra addresses

        std::cout << "Current bra: " << bra.asString() << std::endl;
        std::cout << "I  " << I << std::endl;

        // Annihilate the bra on the bra indices
        if (!bra.annihilateAll(bra_indices, sign)) {  // if we can't annihilate, this doesn't change the bra

            std::cout << "Can't annihilate all bra indices on: " << bra.asString() << std::endl;


            // Go to the beginning of the outer while loop with the next bra
            if (I < dim-1) {  // prevent the last permutation to occur
                fock_space.setNext(bra);
                I++;
                sign = 1;
                std::cout << "I  " << I << std::endl;

                continue;
            } else {
                break;  // we have to jump out if we have looped over the whole bra dimension
            }
        }



        ONV ket = fock_space.get_ONV(0);
        size_t J = 0;
        while (J < dim) {

            std::cout << "Current ket: " << ket.asString() << std::endl;

            // Annihilate the ket on the ket indices
            if (!ket.annihilateAll(ket_indices_reversed, sign)) {  // if we can't annihilate, this doesn't change the ket

                std::cout << "Can't annihilate all ket indices on: " << ket.asString() << std::endl;

                // Go to the beginning of this (the inner) while loop with the next bra
                if (J < dim-1) {  // prevent the last permutation to occur
                    fock_space.setNext(ket);
                    J++;
                    std::cout << "J  " << J << std::endl;

                    sign = 1;
                    continue;
                } else {
                    break;  // we have to jump out if we have looped over the whole ket dimension
                }
            }

//            std::cout << "The fully-annihilated bra is now: " << bra.asString() << std::endl;
//            std::cout << "The fully-annihilated ket is now: " << ket.asString() << std::endl;

            if (bra == ket) {
                std::cout << "I, J: " << I << ',' << J << std::endl;
                std::cout << "sign: " << sign << std::endl;
                std::cout << "c(I), c(J): " << coeff(I) << ',' << coeff(J) << std::endl;
                value += sign * coeff(I) * coeff(J);
            }

            // Reset the previous ket annihilations and move to the next ket
            if (J == dim-1) {  // prevent the last permutation to occur
                break;  // out of the J-loop
            }
            ket.createAll(ket_indices_reversed);
            fock_space.setNext(ket);
            sign = 1;
            J++;
            std::cout << "J  " << J << std::endl;
                

//            std::cout << "I'm stuck in the J loop with J = " << J << " and dim= " << dim << std::endl;


        }  // while J loop

        // Reset the previous bra annihilations and move to the next bra
        if (I == dim-1) {  // prevent the last permutation to occur
            break;  // out of the I-loop
        }
        bra.createAll(bra_indices);
        std::cout << "bra after createAll: " << bra.asString() << std::endl;
        fock_space.setNext(bra);
        sign = 1;
        I++;
        std::cout << "I  " << I << std::endl;


    }  // while I loop


    return value;
}





}  // namespace GQCP
