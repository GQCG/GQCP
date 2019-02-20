// This file is part of GQCG-gqcp.
// 
// Copyright (C) 2017-2019  the GQCG developers
// 
// GQCG-gqcp is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// GQCG-gqcp is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-gqcp.  If not, see <http://www.gnu.org/licenses/>.
// 
#include "FockSpace/FrozenFockSpace.hpp"


namespace GQCP {



/*
 *  CONSTRUCTORS
 */

/**
 *  @param K        the total number of orbitals
 *  @param N        the total number of electrons
 *  @param X        the number of frozen orbitals and electrons
 */
FrozenFockSpace::FrozenFockSpace(size_t K, size_t N, size_t X) :
        BaseFockSpace (K, FockSpace::calculateDimension(K-X, N-X)),
        FockPermutator (N),
        active_fock_space (K-X, N-X),
        X (X)
{}

/**
 *  @param fock_space       (to be frozen) full Fock space
 *  @param X                the number of frozen orbitals and electrons
 */
FrozenFockSpace::FrozenFockSpace(const FockSpace& fock_space, size_t X) :
    FrozenFockSpace(fock_space.get_K(), fock_space.get_N(), X)
{}



/*
 *  PUBLIC OVERRIDEN METHODS
 */

/**
 *  @param representation       a representation of an ONV
 *
 *  @return the next bitstring permutation in the frozen Fock space
 *
 *      Examples (first orbital is frozen):
 *          0101 -> 1001
 *         01101 -> 10011
 */
size_t FrozenFockSpace::ulongNextPermutation(size_t representation) const {
    // generate the permutation from the active space, bitshift left X amount of times to remove the frozen orbital indices before passing it to the active space
    size_t sub_permutation = this->active_fock_space.ulongNextPermutation(representation >> this->X);
    // transform the permutation to the frozen core space, by bitshifting right X amount of times and filling the new 0 bits with 1's (by adding 2^X - 1)
    representation <<= X;
    representation += static_cast<size_t>(pow(2,X)-1);
    return representation;
};

/**
 *  @param representation      a representation of an ONV
 *
 *  @return the address (i.e. the ordering number) of the given ONV
 */
size_t FrozenFockSpace::getAddress(size_t representation) const {
    // transform the representation to the sub space, by bitshifting left X amount of times to remove the frozen orbital indices
    // address of the total ONV in the frozen Fock space is identical to that of the sub ONV in the sub Fock space.
    return this->active_fock_space.getAddress(representation >> this->X);
};

/**
  *  Calculate unsigned representation for a given address
  *
  *  @param address                 the address of the representation is calculated
  *
  *  @return unsigned representation of the address
  */
size_t FrozenFockSpace::calculateRepresentation(size_t address) const {
    // generate the representation in the active space
    size_t representation = this->active_fock_space.calculateRepresentation(address);

    // transform the permutation to the frozen core space, by bitshifting right X amount of times and filling the new 0 bits with 1's (by adding 2^X - 1)
    representation <<= X;
    representation += static_cast<size_t>(pow(2,X)-1);
    return representation;
};



/*
 *  PUBLIC METHODS
 */

/**
 *  @param onv       the ONV
 *
 *  @return the amount of ONVs (with a larger address) this ONV would couple with given a one electron operator
 */
size_t FrozenFockSpace::countOneElectronCouplings(const ONV& onv) const {
    size_t V = K-N;  // amount of virtual orbitals
    size_t coupling_count = 0;

    for (size_t e1 = this->X; e1 < this->N; e1++) {  // start from X to ignore the frozen electrons
        size_t p = onv.get_occupation_index(e1);
        coupling_count += (V + e1 - p);  // amount of virtuals with an index larger than p
    }

    return coupling_count;
}


/**
 *  @param onv       the ONV
 *
 *  @return the amount of ONVs (with a larger address) this ONV would couple with given a two electron operator
 */
size_t FrozenFockSpace::countTwoElectronCouplings(const ONV& onv) const {

    size_t V = K-N; // amount of virtual orbitals
    size_t coupling_count = 0;

    for (size_t e1 = this->X; e1 < this->N; e1++){  // start from X to ignore the frozen electrons

        size_t p = onv.get_occupation_index(e1);
        coupling_count += (V + e1 - p);  //  one electron part

        for (size_t e2 = e1+1; e2 < this->N; e2++){

            size_t q = onv.get_occupation_index(e2);
            size_t coupling_count2 = (V + e2 - q);
            coupling_count += (V-coupling_count2)*coupling_count2;

            if(coupling_count2 > 1 ){
                coupling_count += FockSpace::calculateDimension(coupling_count2, 2);
            }
        }
    }

    return coupling_count;
}


/**
 *  @return the amount non-zero couplings of a one electron coupling scheme in the Fock space
 */
size_t FrozenFockSpace::countTotalOneElectronCouplings() const {
    return this->active_fock_space.countTotalOneElectronCouplings();
}


/**
 *  @return the amount non-zero couplings of a two electron coupling scheme in the Fock space
 */
size_t FrozenFockSpace::countTotalTwoElectronCouplings() const {
    return this->active_fock_space.countTotalTwoElectronCouplings();
}


}  // namespace GQCP
