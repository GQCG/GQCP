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
 *  The full Fock space for a number of orbitals and number of electrons
 *
 *  The ONVs and addresses are linked with a hashing function calculated with an addressing scheme. The implementation of the addressing scheme is from Molecular Electronic-Structure Theory (August 2000) by Trygve Helgaker, Poul Jorgensen, and Jeppe Olsen
 *
 */
class FockSpace: public GQCP::BaseFockSpace {
private:
    const size_t N;  // number of electrons
    Matrixu vertex_weights;  // vertex_weights of the addressing scheme


    // PRIVATE METHODS
    /**
     *  @param representation       a representation of an ONV
     *
     *  @return the next bitstring permutation
     *
     *      Examples:
     *          011 -> 101
     *          101 -> 110
     */
    size_t ulongNextPermutation(size_t representation);


public:
    // CONSTRUCTORS
    /**
     *  @param K        the number of orbitals
     *  @param N        the number of electrons
     */
    FockSpace(size_t K, size_t N);


    // DESTRUCTORS
    ~FockSpace() override = default;


    // GETTERS
    size_t get_vertex_weights(size_t p, size_t m) const { return this->vertex_weights[p][m]; }
    const Matrixu& get_vertex_weights() const { return this->vertex_weights; }
    size_t get_N() const { return this->N; }
    FockSpaceType get_type() const override { return FockSpaceType::FockSpace; }


    // STATIC PUBLIC METHODS
    /**
     *  @param K        the number of orbitals
     *  @param N        the number of electrons
     *
     *  @return the dimension of the Fock space
     */
    static size_t calculateDimension(size_t K, size_t N);


    // PUBLIC METHODS
    /**
     *  @param address      the address (i.e. the ordening number) of the ONV
     *
     *  @return the ONV with the corresponding address
     */
    ONV get_ONV(size_t address);


    /**
     *  Set the current ONV to the next ONV: performs ulongNextPermutation() and updates the corresponding occupation indices of the ONV occupation array
     *
     *  @param onv      the current ONV
     */
    void setNext(ONV& onv);

    /**
     *  @param onv      the ONV
     *
     *  @return the address (i.e. the ordering number) of the given ONV
     */
    size_t getAddress(const ONV& onv);

    /**
     *  Find the next unoccupied orbital in a given ONV,
     *  update the electron count and orbital index,
     *  and calculate a shift in address
     *  resulting from a difference between the initial vertex weights for the encountered occupied orbitals
     *  and the corrected vertex weights accounting for previously annihilated electrons
     *
     *  @param onv       the ONV for which we search the next unnocupied orbital
     *  @param q         the orbital index
     *  @param e         the electron count
     *  @param a         the annihilation count
     *
     *  @return the shift in address resulting from the difference in the corrected electron weights
     */
    size_t shiftUntilNextUnoccupiedOrbital(const ONV& onv, size_t& q, size_t& e, size_t a) const;
};


}  // namespace GQCP


#endif  // GQCP_FOCKSPACE_HPP
