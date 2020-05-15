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


#include "Mathematical/Representation/BlockRankFourTensor.hpp"


namespace GQCP {


template <typename _Scalar>
class T2Amplitudes {
public:
    using Scalar = _Scalar;

private:
    size_t N;                       // the number of occupied orbitals
    size_t M;                       // the total number of orbitals
    BlockRankFourTensor<Scalar> t;  // the T2-amplitudes as a block rank-four tensor, implementing easy operator(i,j,a,b) calls


public:
    /*
     *  CONSTRUCTORS
     */

    /**
     *  Construct the T2-amplitudes given their representation as a BlockRankFourTensor.
     * 
     *  @param t                the T2-amplitudes as a block rank-four tensor, implementing easy operator(i,j,a,b) calls
     *  @param N                the number of occupied orbitals
     *  @param M                the total number of orbitals
     */
    T2Amplitudes(const BlockRankFourTensor<double>& t, const size N, const size_t M) :
        t {t},
        N {N},
        M {M} {}


    /*
     *  OPERATORS
     */

    /**
     *  @param i            an occupied index
     *  @param j            an occupied index
     *  @param a            a virtual index
     *  @param b            a virtual index
     * 
     *  @return the T2-amplitude corresponding to t_{ij}^{ab}
     */
    double operator()(const size_t i, const size_t j, const size_t a, const size_t b) const { return this->t(i, j, a, b); }


    /*
     *  PUBLIC METHODS
     */

    /**
     *  @return the T2-amplitudes as a BlockRankFourTensor
     */
    const BlockRankFourTensor<double>& asBlockRankFourTensor() const { return this->t; }
};


}  // namespace GQCP
