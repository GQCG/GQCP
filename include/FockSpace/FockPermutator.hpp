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
#ifndef GQCP_FOCKPERMUTATOR_HPP
#define GQCP_FOCKPERMUTATOR_HPP


#include "ONV.hpp"


namespace GQCP {


/**
 *  An interface for Fock spaces who manipulate (generate, permutate) ONVs
 */

/**
 *  This class utilizes the CRTP or curiously recurring template pattern.
 *  This allows compile-time identification of pure virtual functions called within a non virtual base method.
 *  This method requires derived classes to derive from the base with the derived class as template parameter
 *  e.g. "Derived : public Base<Derived>"
 *
 *  Example of why this is required is found in this class:
 *       ulongNextPermutation(size_t representation) is pure virtual
 *       setNextONV(ONV &onv) is implemented in the base and calls "ulongNextPermutation(size_t representation)"
 *       and contains "auto& fock_space = static_cast<const DerivedPermutator&>(*this)" a self-cast to the derived instance
 *       allowing compile-time identification of the called method and thus inlining
 */
template<typename DerivedPermutator>
class FockPermutator {
protected:
    size_t N;

public:
    // CONSTRUCTOR
    FockPermutator(size_t N) : N (N)
    {};


    // DESTRUCTOR
    /**
     *  Provide a virtual destructor to make the class abstract
     */
    virtual ~FockPermutator() {};


    // GETTERS
    size_t get_N() const { return this->N; }


    // VIRTUAL PUBLIC METHODS
    /**
     *  @param representation       a representation of an ONV
     *
     *  @return the next bitstring permutation in the Fock space
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


    // PUBLIC METHODS
    /**
     *  @param address      the address (i.e. the ordening number) of the ONV
     *
     *  @return the ONV with the corresponding address
     */
    ONV makeONV(size_t address) const {

        auto& fock_space = static_cast<const DerivedPermutator&>(*this);  // cast to DerivedPermutator for compile-time identification of the function called

        ONV onv (fock_space.get_K(), this->N);
        fock_space.transformONV(onv, address);

        return onv;
    };

    /**
     *  Set the current ONV to the next ONV: performs ulongNextPermutation() and updates the corresponding occupation indices of the ONV occupation array
     *
     *  @param onv      the current ONV
     */
    void setNextONV(ONV &onv) const {
        auto& fock_space = static_cast<const DerivedPermutator&>(*this);  // cast to DerivedPermutator for compile-time identification of the function called
        onv.set_representation(fock_space.ulongNextPermutation(onv.get_unsigned_representation()));
    };

    /**
     *  @param onv      the ONV
     *
     *  @return the address (i.e. the ordering number) of the given ONV
     */
    size_t getAddress(const ONV &onv) const {
        auto& fock_space = static_cast<const DerivedPermutator&>(*this);  // cast to DerivedPermutator for compile-time identification of the function called
        return fock_space.getAddress(onv.get_unsigned_representation());
    };

    /**
     *  Transform an ONV to one corresponding to the given address
     *
     *  @param onv          the ONV
     *  @param address      the address to which the ONV will be set
     */
    void transformONV(ONV &onv, size_t address) const {
        auto& fock_space = static_cast<const DerivedPermutator&>(*this);  // cast to DerivedPermutator for compile-time identification of the function called
        onv.set_representation((fock_space.calculateRepresentation(address)));
    };
};


}  // namespace GQCP


#endif  // GQCP_FOCKPERMUTATOR_HPP
