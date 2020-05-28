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
    SpinUnresolvedFrozenONVBasis(const size_t K, const size_t N, const size_t X);

    /**
     *  @param onv_basis        (to be frozen) full spin-unresolved ONV basis
     *  @param X                the number of frozen orbitals and electrons
     */
    SpinUnresolvedFrozenONVBasis(const SpinUnresolvedONVBasis& onv_basis, const size_t X);


    // DESTRUCTOR

    /**
     *  The default destructor.
     */
    ~SpinUnresolvedFrozenONVBasis() override = default;


    // OVERRIDEN PUBLIC METHODS

    /**
     *  @param representation      a representation of a spin-unresolved ONV
     *
     *  @return the address (i.e. the ordering number) of the given spin-unresolved ONV
     */
    size_t addressOf(const size_t representation) const override;

    /**
     *  @param onv       the spin-unresolved ONV
     *
     *  @return the number of ONVs (with a larger address) this spin-unresolved ONV would couple with given a one electron operator
     */
    size_t countOneElectronCouplings(const SpinUnresolvedONV& onv) const override;

    /**
     *  @return the number of non-zero (non-diagonal) couplings of a one electron coupling scheme in the spin-unresolved ONV basis
     */
    size_t countTotalOneElectronCouplings() const override { return this->active_onv_basis.countTotalOneElectronCouplings(); }

    /**
     *  @return the number of non-zero (non-diagonal) couplings of a two electron coupling scheme in the spin-unresolved ONV basis
     */
    size_t countTotalTwoElectronCouplings() const override { return this->active_onv_basis.countTotalTwoElectronCouplings(); }

    /**
     *  @param onv       the spin-unresolved ONV
     *
     *  @return the number of ONVs (with a larger address) this spin-unresolved ONV would couple with given a two electron operator
     */
    size_t countTwoElectronCouplings(const SpinUnresolvedONV& onv) const override;

    /**
     *  @param representation       a representation of a spin-unresolved ONV
     *
     *  @return the next bitstring permutation in the frozen spin-unresolved ONV basis
     *
     *      Examples (first orbital is frozen):
     *          0101 -> 1001
     *         01101 -> 10011
     */
    size_t nextPermutationOf(const size_t representation) const override;

    /**
      *  Calculate unsigned representation for a given address
      *
      *  @param address                 the address of the representation is calculated
      *
      *  @return unsigned representation of the address
      */
    size_t representationOf(const size_t address) const override;

    /**
     *  @return the type of this ONV basis
     */
    ONVBasisType type() const override { return ONVBasisType::SpinUnresolvedFrozenONVBasis; }


    // PUBLIC METHODS

    /**
     *  @return this' ONV basis that is considered active
     */
    const SpinUnresolvedONVBasis& activeONVBasis() const { return this->active_onv_basis; }

    using ONVManipulator<SpinUnresolvedFrozenONVBasis>::addressOf;

    /**
     *  @return the number of orbitals that are considered frozen in this ONV basis
     */
    size_t numberOfFrozenOrbitals() const { return this->X; }
};


}  // namespace GQCP
