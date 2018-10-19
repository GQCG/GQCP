// This file is part of GQCG-gqcp.
// 
// Copyright (C) 2017-2018  the GQCG developers
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
#ifndef GQCP_FOCKSPACE_HPP
#define GQCP_FOCKSPACE_HPP


#include "FockSpace/BaseFockSpace.hpp"

#include <boost/numeric/conversion/converter.hpp>
#include <boost/math/special_functions.hpp>


namespace GQCP {


/**
 *  The full Fock space for a given set of orbitals and number of electrons
 *  where the ONVs and addresses are linked
 *  through a hashing function calculated with an addressing scheme.
 *  Implementation of the addressing scheme from :
 *      Molecular Electronic-Structure Theory (August 2000) by Trygve Helgaker, Poul Jorgensen, and Jeppe Olsen
 */
class FockSpace: public GQCP::BaseFockSpace {
private:
    const size_t N;  // number of electrons
    Matrixu vertex_weights;  // vertex_weights of the addressing scheme


    // PRIVATE METHODS
    /**
     *  @returns a permutation of the representation, giving the next bitstring permutation in reverse lexical ordering.
     *
     *      Examples:
     *          011 -> 101
     *          101 -> 110
     */
    size_t ulongNextPermutation(size_t representation);


public:
    // CONSTRUCTORS
    /**
     *  Constructor given a @param K (spatial orbitals), N (electrons)
     *  on which the dimensions of the Fock space are based
     */
    FockSpace(size_t K, size_t N);


    // DESTRUCTORS
    ~FockSpace() override = default;


    // GETTERS
    size_t get_vertex_weights(size_t p, size_t m) const { return this->vertex_weights[p][m]; }
    Matrixu get_vertex_weights() const { return this->vertex_weights; }
    size_t get_N() const { return this->N; }
    FockSpaceType get_type() const override { return FockSpaceType::FockSpace; }

    // STATIC PUBLIC METHODS
    /**
     *  Given a number of spatial orbitals @param K
     *  and a number of electrons  @param N,
     *  @return the dimension of the Fock space
     */
    static size_t calculateDimension(size_t K, size_t N);


    // OVERRIDDEN PUBLIC METHODS
    /**
     *  @return the ONV with the corresponding address in the considered space
     */
    ONV get_ONV(size_t address);

    /**
     *  sets @param ONV to the next ONV in the space
     *  performs the ulongNextPermutation() function
     *  and updates the corresponding occupation indexes
     *  of the ONV occupation array
     */
    void setNext(ONV& onv);

    /**
     *  @return the Fock space address (i.e. the ordering number) of the @param onv in reverse lexical ordering, in the fock space.
     */
    size_t getAddress(const ONV& onv);


    // FRIEND CLASSES
    friend class DOCI;
    friend class FCI;
};


}  // namespace GQCP


#endif  // GQCP_FOCKSPACE_HPP
