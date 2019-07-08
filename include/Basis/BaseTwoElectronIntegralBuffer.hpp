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
 *  @tparam _Scalar         the scalar representation of an integral
 *  @tparam _N              the number of components the operator has
 */
template <typename _Scalar, size_t _N>
class BaseTwoElectronIntegralBuffer {
public:
    using Scalar = _Scalar;  // the scalar representation of an integral
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
     *  @return the matrix representation of the integrals that are in this buffer
     */
    virtual std::array<Tensor<Scalar, 4>, N> integrals() const = 0;



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
     *  Place the calculated integrals over the shells inside the total matrix representation
     * 
     *  @param full_components          the components of the full matrix representation (over all the basis functions) of the operator
     *  @param bf1                      the total basis function index of the first basis function in the first shell
     *  @param bf2                      the total basis function index of the first basis function in the second shell
     *  @param bf3                      the total basis function index of the first basis function in the third shell
     *  @param bf4                      the total basis function index of the first basis function in the fourth shell     
     */
    void emplace(std::array<SquareRankFourTensor<Scalar>, N>& full_components, const size_t bf1, const size_t bf2, const size_t bf3, const size_t bf4) const {

        auto components = this->integrals();  // N components
        for (size_t i = 0; i < N; i++) {
            auto& full_component = full_components[i];
            auto& component = components[i];

            // Place the calculated integrals inside the correct block (for Tensors, this is a 'slice')
            const auto bf1_int = static_cast<int>(bf1);
            const auto bf2_int = static_cast<int>(bf2);
            const auto bf3_int = static_cast<int>(bf3);
            const auto bf4_int = static_cast<int>(bf4);

            const auto nbf1_int = static_cast<int>(this->nbf1);
            const auto nbf2_int = static_cast<int>(this->nbf2);
            const auto nbf3_int = static_cast<int>(this->nbf3);
            const auto nbf4_int = static_cast<int>(this->nbf4);

            Eigen::array<int, 4> offsets {bf1_int, bf2_int, bf3_int, bf4_int};
            Eigen::array<int, 4> extents {nbf1_int, nbf2_int, nbf3_int, nbf4_int};  // length of the slices
            full_component.Eigen().slice(offsets, extents) = component.Eigen();
        }
    }
};


}  // namespace GQCP


#endif  // GQCP_BASETWOELECTRONINTEGRALBUFFER_HPP
