#include "RDM/NRDMCalculator.hpp"



#include <iostream>



namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param fock_space       the FCI Fock space with spin orbitals
 *  @param coeff            the expansion coefficient vector
 */
NRDMCalculator::NRDMCalculator(const FockSpace& fock_space, const Eigen::VectorXd& coeff) :
    fock_space (fock_space),
    coeff (coeff)
{}



/*
 *  PUBLIC METHODS
 */

/**
 *  @param bra_indices      the indices of the orbitals that should be annihilated on the left (on the bra)
 *  @param ket_indices      the indices of the orbitals that should be annihilated on the right (on the ket)
 *
 *  @return an element of the N-RDM, as specified by the given bra and ket indices
 *
 *      calculateElement({0, 1}, {2, 1}) would calculate d^{(2)} (0, 1, 1, 2)
 */
double NRDMCalculator::calculateElement(const std::vector<size_t>& bra_indices, const std::vector<size_t>& ket_indices) const {


    // The ket indices should be reversed because the annihilators on the ket should be applied from right to left
    std::vector<size_t> ket_indices_reversed = ket_indices;
    std::reverse(ket_indices_reversed.begin(), ket_indices_reversed.end());


    double value = 0.0;
    int sign = 1;
    FockSpace fock_space = this->fock_space;  // make a copy because this method is marked const
    size_t dim = fock_space.get_dimension();


    ONV bra = fock_space.get_ONV(0);
    size_t I = 0;
    while (I < dim) {  // loop over all bra addresses

        // Annihilate the bra on the bra indices
        if (!bra.annihilateAll(bra_indices, sign)) {  // if we can't annihilate, the bra doesn't change

            // Go to the beginning of the outer while loop with the next bra
            if (I < dim-1) {  // prevent the last permutation to occur
                fock_space.setNext(bra);
                I++;
                sign = 1;
                continue;
            } else {
                break;  // we have to jump out if we have looped over the whole bra dimension
            }
        }


        ONV ket = fock_space.get_ONV(0);
        size_t J = 0;
        while (J < dim) {  // loop over all ket indices

            // Annihilate the ket on the ket indices
            if (!ket.annihilateAll(ket_indices_reversed, sign)) {  // if we can't annihilate, the ket doesn't change

                // Go to the beginning of this (the inner) while loop with the next bra
                if (J < dim-1) {  // prevent the last permutation to occur
                    fock_space.setNext(ket);
                    J++;
                    sign = 1;
                    continue;
                } else {
                    break;  // we have to jump out if we have looped over the whole ket dimension
                }
            }

            if (bra == ket) {
                value += sign * this->coeff(I) * this->coeff(J);
            }

            // Reset the previous ket annihilations and move to the next ket
            if (J == dim-1) {  // prevent the last permutation to occur
                break;  // out of the J-loop
            }
            ket.createAll(ket_indices_reversed);
            fock_space.setNext(ket);
            sign = 1;
            J++;
        }  // while J loop

        // Reset the previous bra annihilations and move to the next bra
        if (I == dim-1) {  // prevent the last permutation to occur
            break;  // out of the I-loop
        }
        bra.createAll(bra_indices);
        fock_space.setNext(bra);
        sign = 1;
        I++;
    }  // while I loop

    return value;
}



}  // namespace GQCP
