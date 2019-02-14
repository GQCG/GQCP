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
 *  @param K        the number of orbitals
 *  @param N        the number of electrons
 *  @param X        the number of frozen orbitals
 */
FrozenFockSpace::FrozenFockSpace(size_t K, size_t N, size_t X) :
        BaseFockSpace(K, FockSpace::calculateDimension(K-X, N-X)),
        fock_space(K-X, N-X),
        N (N),
        X (X)
{}

/**
 *  @param fock_space       non-frozen sub Fock space
 *  @param X                the number of frozen orbitals
 */
FrozenFockSpace::FrozenFockSpace(const FockSpace& fock_space, size_t X) :
    FrozenFockSpace(fock_space.get_K(), fock_space.get_N(), X)
{}



/*
 *  PUBLIC METHODS
 */

/**
 *  @param address      the address (i.e. the ordening number) of the ONV
 *
 *  @return the ONV with the corresponding address
 */
ONV FrozenFockSpace::makeONV(size_t address) const {
    ONV onv (this->K, this->N);
    this->transformONV(onv, address);
    return onv;
}


/**
 *  Set the current ONV to the next ONV
 *
 *  @param onv      the current ONV
 */
void FrozenFockSpace::setNextONV(ONV& onv) const {
    onv.set_representation((this->fock_space.ulongNextPermutation(onv.get_unsigned_representation() >> X) << X) + pow(2,X)-1);
}


/**
 *  @param onv                                              the ONV
 *
 *  @return the address (i.e. the ordering number) of the given ONV
 */
size_t FrozenFockSpace::getAddress(const ONV& onv) const {
    return this->fock_space.getAddress(onv.get_unsigned_representation() >> X);
}


/**
 *  Transform an ONV to one with corresponding to the given address
 *
 *  @param onv          the ONV
 *  @param address      the address to which the ONV will be set
 */
void FrozenFockSpace::transformONV(ONV& onv, size_t address) const {
    onv.set_representation((this->fock_space.calculateRepresentation(address) << X) + pow(2,X)-1);
}


/**
 *  @param onv       the ONV
 *
 *  @return the amount of ONVs (with a larger address) this ONV would couple with given a one electron operator
 */
size_t FrozenFockSpace::countOneElectronCouplings(const ONV& onv) const {
    size_t V = K-N;  // amount of virtual orbitals
    size_t coupling_count = 0;

    for (size_t e1 = 1; e1 < this->N; e1++) {
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

    for (size_t e1 = 1; e1 < this->N; e1++){

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
    return this->fock_space.countTotalOneElectronCouplings();
}


/**
 *  @return the amount non-zero couplings of a two electron coupling scheme in the Fock space
 */
size_t FrozenFockSpace::countTotalTwoElectronCouplings() const {
    return this->fock_space.countTotalTwoElectronCouplings();
}


}  // namespace GQCP
