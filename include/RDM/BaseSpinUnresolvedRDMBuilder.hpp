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
#ifndef GQCG_BASESPINUNRESOLVEDRDMBUILDER_HPP
#define GQCG_BASESPINUNRESOLVEDRDMBUILDER_HPP


#include "RDM/OneRDM.hpp"
#include "RDM/TwoRDM.hpp"
#include "RDM/RDMs.hpp"
#include "FockSpace/BaseFockSpace.hpp"


namespace GQCP {


/**
 *  BaseSpinUnresolvedRDMBuilder is an abstract base class whose derived classes are capable of calculating 1- and 2-RDM from wave functions in a spin unresolved basis
 */
class BaseSpinUnresolvedRDMBuilder {
public:
    // CONSTRUCTOR
    BaseSpinUnresolvedRDMBuilder() = default;


    // DESTRUCTOR
    /**
     *  Provide a pure virtual destructor to make the class abstract
     */
    virtual ~BaseSpinUnresolvedRDMBuilder() = 0;


    // PURE VIRTUAL GETTERS
    virtual BaseFockSpace* get_fock_space() = 0;


    // PURE VIRTUAL PUBLIC METHODS
    /**
     *  @param x        the coefficient vector representing the wave function
     *
     *  @return the 1-RDM given a coefficient vector
     */
    virtual OneRDM calculate1RDM(const Eigen::VectorXd& x) const = 0;

    /**
     *  @param x        the coefficient vector representing the wave function
     *
     *  @return the 2-RDM given a coefficient vector
     */
    virtual TwoRDM calculate2RDM(const Eigen::VectorXd& x) const = 0;

    /**
     *  @param bra_indices      the indices of the orbitals that should be annihilated on the left (on the bra)
     *  @param ket_indices      the indices of the orbitals that should be annihilated on the right (on the ket)
     *  @param x                the coefficient vector representing the wave function
     *
     *  @return an element of the N-RDM, as specified by the given bra and ket indices.
     *
     *      calculateElement({0, 1}, {2, 1}) would calculate d^{(2)} (0, 1, 1, 2): the operator string would be a^\dagger_0 a^\dagger_1 a_2 a_1
     */
    virtual double calculateElement(const std::vector<size_t>& bra_indices, const std::vector<size_t>& ket_indices, const Eigen::VectorXd& x) const = 0;
};


}  // namespace GQCG


#endif  // GQCG_BASESPINUNRESOLVEDRDMBUILDER_HPP
