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
#ifndef GQCP_BASETWOELECTRONINTEGRALBUFFER_HPP
#define GQCP_BASETWOELECTRONINTEGRALBUFFER_HPP


#include "Mathematical/SquareRankFourTensor.hpp"

#include <array>


namespace GQCP {


/**
 *  A base class to implement buffers for storing two-electron integrals
 * 
 *  An integral buffer facilitates placing integrals (see emplace()) that were calculated over shells into the correct matrix representation of an operator in an orbital basis
 * 
 *  @tparam _IntegralScalar         the scalar representation of an integral
 *  @tparam _N                      the number of components the operator has
 */
template <typename _IntegralScalar, size_t _N>
class BaseTwoElectronIntegralBuffer {
public:
    using IntegralScalar = _IntegralScalar;  // the scalar representation of an integral
    static constexpr auto N = _N;  // the number of components the operator has


protected:
    size_t nbf1;  // the number of basis functions in the first shell
    size_t nbf2;  // the number of basis functions in the second shell
    size_t nbf3;  // the number of basis functions in the third shell
    size_t nbf4;  // the number of basis functions in the fourth shell



public:

    /*
     *  CONSTRUCTORS
     */

    /**
     *  @param nbf1             the number of basis functions in the first shell
     *  @param nbf2             the number of basis functions in the second shell
     *  @param nbf3             the number of basis functions in the third shell
     *  @param nbf4             the number of basis functions in the fourth shell
     */
    BaseTwoElectronIntegralBuffer(const size_t nbf1, const size_t nbf2, const size_t nbf3, const size_t nbf4) :
        nbf1 (nbf1),
        nbf2 (nbf2),
        nbf3 (nbf3),
        nbf4 (nbf4)
    {}



    /*
     *  PUBLIC PURE VIRTUAL METHODS
     */

    /**
     *  @param i            the index of the component of the operator
     *  @param f1           the index of the basis function within shell 1
     *  @param f2           the index of the basis function within shell 2
     *  @param f3           the index of the basis function within shell 3
     *  @param f4           the index of the basis function within shell 4
     * 
     *  @return a value from this integral buffer
     */
    virtual IntegralScalar value(const size_t i, const size_t f1, const size_t f2, const size_t f3, const size_t f4) const = 0;


    /*
     *  PUBLIC METHODS
     */

    /**
     *  @return the number of basis functions that are in the first shell
     */
    size_t numberOfBasisFunctionsInShell1() const { return this->nbf1; }

    /**
     *  @return the number of basis functions that are in the second shell
     */
    size_t numberOfBasisFunctionsInShell2() const { return this->nbf2; }

    /**
     *  @return the number of basis functions that are in the third shell
     */
    size_t numberOfBasisFunctionsInShell3() const { return this->nbf3; }

    /**
     *  @return the number of basis functions that are in the fourth shell
     */
    size_t numberOfBasisFunctionsInShell4() const { return this->nbf4; }

    /**
     *  Place the calculated integrals inside the matrix representation of the integrals
     * 
     *  @param full_components          the components of the full matrix representation (over all the basis functions) of the operator
     *  @param bf1                      the total basis function index of the first basis function in the first shell
     *  @param bf2                      the total basis function index of the first basis function in the second shell
     *  @param bf3                      the total basis function index of the first basis function in the third shell
     *  @param bf4                      the total basis function index of the first basis function in the fourth shell     
     */
    void emplace(std::array<SquareRankFourTensor<IntegralScalar>, N>& full_components, const size_t bf1, const size_t bf2, const size_t bf3, const size_t bf4) const {

        // Place the calculated integrals inside the matrix representation of the integrals
        for (size_t f1 = 0; f1 != this->nbf1; f1++) {  // f1: index of basis function within shell 1
            for (size_t f2 = 0; f2 != this->nbf2; f2++) {  // f2: index of basis function within shell 2
                for (size_t f3 = 0; f3 != this->nbf3; f3++) {  // f3: index of basis function within shell 3
                    for (size_t f4 = 0; f4 != this->nbf4; f4++) {  // f4: index of basis function within shell 4

                        for (size_t i = 0; i < N; i++) {
                            full_components[i](bf1 + f1, bf2 + f2, bf3 + f3, bf4 + f4) = this->value(i, f1, f2, f3, f4);  // in chemist's notation
                        }  // components

                    }  // f1
                }  // f2
            }  // f3
        }  // f4
    }
};


}  // namespace GQCP


#endif  // GQCP_BASETWOELECTRONINTEGRALBUFFER_HPP
