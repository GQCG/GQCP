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
#ifndef GQCP_BASEONEELECTRONINTEGRALBUFFER_HPP
#define GQCP_BASEONEELECTRONINTEGRALBUFFER_HPP


#include "Mathematical/SquareMatrix.hpp"

#include <array>


namespace GQCP {


/**
 *  A base class to implement buffers for storing one-electron integrals
 * 
 *  An integral buffer facilitates placing integrals (see emplace()) that were calculated over shells into the correct matrix representation of an operator in an orbital basis
 * 
 *  @tparam _Scalar         the scalar representation of an integral
 *  @tparam _N              the number of components the operator has
 */
template <typename _Scalar, size_t _N>
class BaseOneElectronIntegralBuffer {
public:
    using Scalar = _Scalar;  // the scalar representation of an integral
    constexpr static auto N = _N;  // the number of components the operator has


private:
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
     *  PUBLIC PURE VIRTUAL METHODS
     */

    /**
     *  @return the matrix representation of the integrals that are in this buffer
     */
    virtual std::array<SquareMatrix<Scalar>, N> matrix() const = 0;



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
     *  Place the calculated integrals over some shells inside the total matrix representation
     * 
     *  @param full_components          the components of the full matrix representation (over all the basis functions) of the operator
     *  @param bf1                      the total basis function index of the first basis function in the first shell
     *  @param bf2                      the total basis function index of the first basis function in the second shell
     */
    void emplace(std::array<SquareMatrix<Scalar>, N>& full_components, const size_t bf1, const size_t bf2) const {

        const auto components = this->matrix();  // N components
        for (size_t i = 0; i < N; i++) {
            auto& full_component = full_components[i];
            const auto& component = components[i];

            full_component.block(bf1, bf2, nbf1, nbf2) = component;
        }
    }
};


}  // namespace GQCP


#endif  // GQCP_BASEONEELECTRONINTEGRALBUFFER_HPP
