// This file is part of GQCG-GQCP.
//
// Copyright (C) 2017-2020  the GQCG developers
//
// GQCG-GQCP is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// GQCG-GQCP is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-GQCP.  If not, see <http://www.gnu.org/licenses/>.

#pragma once


#include "ONVBasis/BaseFrozenCoreONVBasis.hpp"
#include "ONVBasis/BaseONVBasis.hpp"
#include "ONVBasis/ONVManipulator.hpp"
#include "ONVBasis/SpinUnresolvedONVBasis.hpp"


namespace GQCP {


/**
 *  A spin-unresolved frozen ONV basis: this is a subset of the N-electron spin-unresolved ONV basis in which the first X orbitals are always occupied.
 */
class SpinUnresolvedFrozenONVBasis: public BaseFrozenCoreONVBasis, public ONVManipulator<SpinUnresolvedFrozenONVBasis> {
protected:
    size_t X;                                 // number of frozen orbitals/electrons
    SpinUnresolvedONVBasis active_onv_basis;  // active (non-frozen) spin-unresolved ONV basis containing only the active electrons (N-X) and orbitals (K-X)


public:
    // CONSTRUCTORS
    /**
     *  @param K        the total number of orbitals
     *  @param N        the total number of electrons
     *  @param X        the number of frozen orbitals and electrons
     */
    SpinUnresolvedFrozenONVBasis(size_t K, size_t N, size_t X);

    /**
     *  @param onv_basis        (to be frozen) full spin-resolved ONV basis
     *  @param X                the number of frozen orbitals and electrons
     */
    SpinUnresolvedFrozenONVBasis(const SpinUnresolvedONVBasis& onv_basis, size_t X);


    // DESTRUCTOR
    ~SpinUnresolvedFrozenONVBasis() override = default;


    // GETTERS
    size_t get_number_of_frozen_orbitals() const { return this->X; }
    const SpinUnresolvedONVBasis& get_active_fock_space() const { return this->active_onv_basis; }
    ONVBasisType get_type() const override { return ONVBasisType::SpinUnresolvedFrozenONVBasis; }


    // OVERRIDEN PUBLIC METHODS
    /**
     *  @param representation       a representation of a spin-resolved ONV
     *
     *  @return the next bitstring permutation in the frozen spin-resolved ONV basis
     *
     *      Examples (first orbital is frozen):
     *          0101 -> 1001
     *         01101 -> 10011
     */
    size_t ulongNextPermutation(size_t representation) const override;

    /**
     *  @param representation      a representation of a spin-resolved ONV
     *
     *  @return the address (i.e. the ordering number) of the given spin-resolved ONV
     */
    size_t getAddress(size_t representation) const override;

    /**
      *  Calculate unsigned representation for a given address
      *
      *  @param address                 the address of the representation is calculated
      *
      *  @return unsigned representation of the address
      */
    size_t calculateRepresentation(size_t address) const override;

    /**
     *  @param onv       the spin-resolved ONV
     *
     *  @return the number of ONVs (with a larger address) this spin-resolved ONV would couple with given a one electron operator
     */
    size_t countOneElectronCouplings(const SpinUnresolvedONV& onv) const override;

    /**
     *  @param onv       the spin-resolved ONV
     *
     *  @return the number of ONVs (with a larger address) this spin-resolved ONV would couple with given a two electron operator
     */
    size_t countTwoElectronCouplings(const SpinUnresolvedONV& onv) const override;

    /**
     *  @return the amount non-zero (non-diagonal) couplings of a one electron coupling scheme in the spin-resolved ONV basis
     */
    size_t countTotalOneElectronCouplings() const override;

    /**
     *  @return the amount non-zero (non-diagonal) couplings of a two electron coupling scheme in the spin-resolved ONV basis
     */
    size_t countTotalTwoElectronCouplings() const override;


    // PUBLIC METHODS
    /**
     *  If we have
     *      SpinUnresolvedFrozenONVBasis onv_basis;
     *
     *  This makes sure that we can call
     *      onv_basis.getAddress(onv);
     *  instead of the syntax
     *      onv_basis.ONVManipulator<SpinUnresolvedFrozenONVBasis>::getAddress(onv);
     */
    using ONVManipulator<SpinUnresolvedFrozenONVBasis>::getAddress;
};


}  // namespace GQCP
