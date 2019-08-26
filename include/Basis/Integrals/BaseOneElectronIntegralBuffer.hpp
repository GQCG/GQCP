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


#include "Mathematical/ChemicalMatrix.hpp"

#include <array>


namespace GQCP {


/**
 *  A base class to implement buffers for storing one-electron integrals
 * 
 *  An integral buffer facilitates placing integrals (see emplace()) that were calculated over shells into the correct matrix representation of an operator in an orbital basis
 * 
 *  @tparam _IntegralScalar         the scalar representation of an integral
 *  @tparam _N                      the number of components the corresponding operator has
 */
template <typename _IntegralScalar, size_t _N>
class BaseOneElectronIntegralBuffer {
public:
    using IntegralScalar = _IntegralScalar;  // the scalar representation of an integral
    static constexpr auto N = _N;  // the number of components the operator has


protected:
    size_t nbf1;  // the number of basis functions in the first shell
    size_t nbf2;  // the number of basis functions in the second shell


public:

    /*
     *  CONSTRUCTORS
     */

    /**
     *  @param nbf1             the number of basis functions in the first shell
     *  @param nbf2             the number of basis functions in the second shell
     */
    BaseOneElectronIntegralBuffer(const size_t nbf1, const size_t nbf2) :
        nbf1 (nbf1),
        nbf2 (nbf2)
    {}


    /*
     *  DESTRUCTOR
     */
    virtual ~BaseOneElectronIntegralBuffer() = default;



    /*
     *  PUBLIC PURE VIRTUAL METHODS
     */

    /**
     *  @param i            the index of the component of the operator
     *  @param f1           the index of the basis function within shell 1
     *  @param f2           the index of the basis function within shell 2
     * 
     *  @return a value from this integral buffer
     */
    virtual IntegralScalar value(const size_t i, const size_t f1, const size_t f2) const = 0;

    /**
     *  @return if all the values of the calculated integrals are zero
     */
    virtual bool areIntegralsAllZero() const = 0;


    /*
     *  PUBLIC METHODS
     */

    /**
     *  @return the number of basis functions that are in the first shell
     */
    size_t numberOfBasisFunctionsInShell1() const { return this->nbf1; }

    /**
     *  @return the number of basis functions that are in the first shell
     */
    size_t numberOfBasisFunctionsInShell2() const { return this->nbf2; }

    /**
     *  Place the calculated integrals inside the matrix representation of the integrals
     * 
     *  @param full_components          the components of the full matrix representation (over all the basis functions) of the operator
     *  @param bf1                      the total basis function index of the first basis function in the first shell
     *  @param bf2                      the total basis function index of the first basis function in the second shell
     */
    void emplace(std::array<ChemicalMatrix<IntegralScalar>, N>& full_components, const size_t bf1, const size_t bf2) const {

        // Place the calculated integrals inside the matrix representation of the integrals
        for (size_t f1 = 0; f1 < this->nbf1; f1++) {  // the index of the basis function within shell 1
            for (size_t f2 = 0; f2 < this->nbf2; f2++) {  // the index of the basis function within shell 2
                
                for (size_t i = 0; i < N; i++) {  // loop over the components
                    full_components[i](bf1 + f1, bf2 + f2) = this->value(i, f1, f2);
                }  // components

            }  // f2
        }  // f1
    }
};


}  // namespace GQCP
