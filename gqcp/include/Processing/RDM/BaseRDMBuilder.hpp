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


#include "ONVBasis/BaseONVBasis.hpp"
#include "Processing/RDM/OneRDM.hpp"
#include "Processing/RDM/RDMs.hpp"
#include "Processing/RDM/TwoRDM.hpp"


namespace GQCP {


/**
 *  BaseRDMBuilder is an abstract base class whose derived classes are capable of calculating 1- and 2-RDMs from wave functions in a spin resolved basis
 */
class BaseRDMBuilder {
public:
    // CONSTRUCTORS
    BaseRDMBuilder() = default;


    // DESTRUCTOR
    virtual ~BaseRDMBuilder() = default;


    // PURE VIRTUAL GETTERS
    virtual const BaseONVBasis* get_fock_space() const = 0;


    // PURE VIRTUAL PUBLIC METHODS

    /**
     *  @param x        the coefficient vector representing the wave function
     *
     *  @return all 1-RDMs given a coefficient vector
     */
    virtual OneRDMs<double> calculate1RDMs(const VectorX<double>& x) const = 0;

    /**
     *  @param x        the coefficient vector representing the wave function
     *
     *  @return all 2-RDMs given a coefficient vector
     */
    virtual TwoRDMs<double> calculate2RDMs(const VectorX<double>& x) const = 0;

    /**
     *  @param bra_indices      the indices of the orbitals that should be annihilated on the left (on the bra)
     *  @param ket_indices      the indices of the orbitals that should be annihilated on the right (on the ket)
     *  @param x                the coefficient vector representing the wave function
     *
     *  @return an element of the spin-summed (total) N-RDM, as specified by the given bra and ket indices.
     *
     *      calculateElement({0, 1}, {2, 1}) would calculate d^{(2)} (0, 1, 1, 2): the operator string would be a^\dagger_0 a^\dagger_1 a_2 a_1
     */
    virtual double calculateElement(const std::vector<size_t>& bra_indices, const std::vector<size_t>& ket_indices, const VectorX<double>& x) const = 0;
};


}  // namespace GQCP
