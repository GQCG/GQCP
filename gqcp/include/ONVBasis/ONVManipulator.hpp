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
#pragma once


#include "ONVBasis/ONV.hpp"
#include "Utilities/CRTP.hpp"


namespace GQCP {


/**
 *  An interface for ONVBases that generate, permutate and can calculate the address of an ONV given its unsigned representation.
 * 
 *  @tparam _DerivedManipulator             the type of the derived class
 */
template <typename _DerivedManipulator>
class ONVManipulator : public CRTP<_DerivedManipulator> {
protected:
    size_t N;  // number of electrons

public:

    /*
     *  CONSTRUCTORS
     */

    /**
     *  Default constructor setting everything to zero
     */
    ONVManipulator() :
        ONVManipulator(0)
    {}

    /**
     *  @param N            the number of electrons
     */
    ONVManipulator(const size_t N) :
        N (N)
    {}


    // DESTRUCTOR

    virtual ~ONVManipulator() = default;


    // GETTERS
    size_t get_N() const { return this->N; }


    // VIRTUAL PUBLIC METHODS

    /**
     *  @param representation       a representation of an ONV
     *
     *  @return the next bitstring permutation in the ONV basis
     */
    virtual size_t ulongNextPermutation(size_t representation) const = 0;

    /**
     *  @param representation      a representation of an ONV
     *
     *  @return the address (i.e. the ordering number) of the given ONV
     */
    virtual size_t getAddress(size_t representation) const = 0;

    /**
      *  Calculate unsigned representation for a given address
      *
      *  @param address                 the address of the representation is calculated
      *
      *  @return unsigned representation of the address
      */
    virtual size_t calculateRepresentation(size_t address) const = 0;

    /**
     *  @param onv       the ONV
     *
     *  @return the amount of ONVs (with a larger address) this ONV would couple with given a one electron operator
     */
    virtual size_t countOneElectronCouplings(const ONV& onv) const = 0;

    /**
     *  @param onv       the ONV
     *
     *  @return the amount of ONVs (with a larger address) this ONV would couple with given a two electron operator
     */
    virtual size_t countTwoElectronCouplings(const ONV& onv) const = 0;

    /**
     *  @return the amount non-zero (non-diagonal) couplings of a one electron coupling scheme in the ONV basis
     */
    virtual size_t countTotalOneElectronCouplings() const = 0;

    /**
     *  @return the amount non-zero (non-diagonal) couplings of a two electron coupling scheme in the ONV basis
     */
    virtual size_t countTotalTwoElectronCouplings() const = 0;


    // PUBLIC METHODS

    /**
     *  @param address      the address (i.e. the ordening number) of the ONV
     *
     *  @return the ONV with the corresponding address
     */
    ONV makeONV(size_t address) const {

        const auto& fock_space = this->derived();

        ONV onv (fock_space.get_K(), this->N);
        fock_space.transformONV(onv, address);
        return onv;
    }

    /**
     *  Set the current ONV to the next ONV: performs ulongNextPermutation() and updates the corresponding occupation indices of the ONV occupation array
     *
     *  @param onv      the current ONV
     */
    void setNextONV(ONV &onv) const {

        const auto& fock_space = this->derived();
        onv.set_representation(fock_space.ulongNextPermutation(onv.get_unsigned_representation()));
    }

    /**
     *  @param onv      the ONV
     *
     *  @return the address (i.e. the ordering number) of the given ONV
     */
    size_t getAddress(const ONV &onv) const {

        const auto& fock_space = this->derived();
        return fock_space.getAddress(onv.get_unsigned_representation());
    };

    /**
     *  Transform an ONV to one corresponding to the given address
     *
     *  @param onv          the ONV
     *  @param address      the address to which the ONV will be set
     */
    void transformONV(ONV &onv, size_t address) const {

        const auto& fock_space = this->derived();
        onv.set_representation((fock_space.calculateRepresentation(address)));
    }
};


}  // namespace GQCP
