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
#include "ONVBasis/SpinUnresolvedFrozenONVBasis.hpp"


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param K        the total number of orbitals
 *  @param N        the total number of electrons
 *  @param X        the number of frozen orbitals and electrons
 */
SpinUnresolvedFrozenONVBasis::SpinUnresolvedFrozenONVBasis(size_t K, size_t N, size_t X) :
    BaseFrozenCoreONVBasis(std::make_shared<SpinUnresolvedONVBasis>(SpinUnresolvedONVBasis(K - X, N - X)), X),
    ONVManipulator(N),
    active_onv_basis {K - X, N - X},
    X {X} {}


/**
 *  @param fock_space       (to be frozen) full spin-unresolved basis
 *  @param X                the number of frozen orbitals and electrons
 */
SpinUnresolvedFrozenONVBasis::SpinUnresolvedFrozenONVBasis(const SpinUnresolvedONVBasis& fock_space, size_t X) :
    SpinUnresolvedFrozenONVBasis(fock_space.get_K(), fock_space.get_N(), X) {}


/*
 *  PUBLIC OVERRIDEN METHODS
 */

/**
 *  @param representation       a representation of an ONV
 *
 *  @return the next bitstring permutation in this ONV basis
 *
 *      Examples (first orbital is frozen):
 *          0101 -> 1001
 *         01101 -> 10011
 */
size_t SpinUnresolvedFrozenONVBasis::ulongNextPermutation(size_t representation) const {

    // Generate the permutation from the active space, bitshift left X number of times to remove the frozen orbital indices before passing it to the active space
    size_t sub_permutation = this->active_onv_basis.ulongNextPermutation(representation >> this->X);

    // Transform the permutation to the frozen core space, by bitshifting right X number of times and filling the new 0 bits with 1's (by adding 2^X - 1)
    sub_permutation <<= X;
    sub_permutation += static_cast<size_t>(pow(2, X) - 1);
    return sub_permutation;
};


/**
 *  @param representation      a representation of an ONV
 *
 *  @return the address (i.e. the ordering number) of the given ONV
 */
size_t SpinUnresolvedFrozenONVBasis::getAddress(size_t representation) const {

    // Transform the representation to the sub space, by bitshifting left X number of times to remove the frozen orbital indices
    // Address of the total ONV in this frozen ONV basis is identical to that of the sub ONV in the sub ONV basis.
    return this->active_onv_basis.getAddress(representation >> this->X);
};


/**
  *  Calculate unsigned representation for a given address
  *
  *  @param address                 the address of the representation is calculated
  *
  *  @return unsigned representation of the address
  */
size_t SpinUnresolvedFrozenONVBasis::calculateRepresentation(size_t address) const {
    // generate the representation in the active space
    size_t representation = this->active_onv_basis.calculateRepresentation(address);

    // transform the permutation to the frozen core space, by bitshifting right X number of times and filling the new 0 bits with 1's (by adding 2^X - 1)
    representation <<= X;
    representation += static_cast<size_t>(pow(2, X) - 1);
    return representation;
};


/*
 *  PUBLIC METHODS
 */

/**
 *  @param onv       the ONV
 *
 *  @return the number of ONVs (with a larger address) the given ONV would couple with given a one electron operator
 */
size_t SpinUnresolvedFrozenONVBasis::countOneElectronCouplings(const SpinUnresolvedONV& onv) const {

    size_t V = K - N;  // number of virtual orbitals
    size_t coupling_count = 0;

    for (size_t e1 = this->X; e1 < this->N; e1++) {  // start from X to ignore the frozen electrons
        size_t p = onv.get_occupation_index(e1);
        coupling_count += (V + e1 - p);  // number of virtuals with an index larger than p
    }

    return coupling_count;
}


/**
 *  @param onv       the ONV
 *
 *  @return the number of ONVs (with a larger address) the given ONV would couple with given a two electron operator
 */
size_t SpinUnresolvedFrozenONVBasis::countTwoElectronCouplings(const SpinUnresolvedONV& onv) const {

    size_t V = K - N;  // number of virtual orbitals
    size_t coupling_count = 0;

    for (size_t e1 = this->X; e1 < this->N; e1++) {  // start from X to ignore the frozen electrons

        size_t p = onv.get_occupation_index(e1);
        coupling_count += (V + e1 - p);  //  one electron part

        for (size_t e2 = e1 + 1; e2 < this->N; e2++) {

            size_t q = onv.get_occupation_index(e2);
            size_t coupling_count2 = (V + e2 - q);
            coupling_count += (V - coupling_count2) * coupling_count2;

            if (coupling_count2 > 1) {
                coupling_count += SpinUnresolvedONVBasis::calculateDimension(coupling_count2, 2);
            }
        }
    }

    return coupling_count;
}


/**
 *  @return the amount non-zero couplings of a one electron coupling scheme in this ONV basis
 */
size_t SpinUnresolvedFrozenONVBasis::countTotalOneElectronCouplings() const {
    return this->active_onv_basis.countTotalOneElectronCouplings();
}


/**
 *  @return the amount non-zero couplings of a two electron coupling scheme in this ONV basis
 */
size_t SpinUnresolvedFrozenONVBasis::countTotalTwoElectronCouplings() const {
    return this->active_onv_basis.countTotalTwoElectronCouplings();
}


}  // namespace GQCP
