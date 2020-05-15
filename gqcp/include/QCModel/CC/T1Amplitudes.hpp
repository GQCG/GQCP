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


#include "Mathematical/Representation/BlockMatrix.hpp"


namespace GQCP {


/**
 *  The coupled-cluster T1 amplitudes. According to context, this class may either represent restricted (i.e. spatial-orbital) amplitudes, or generalized (spinor) amplitudes.
 */
template <typename _Scalar>
class T1Amplitudes {
public:
    using Scalar = _Scalar;


private:
    size_t N;               // the number of occupied orbitals
    size_t M;               // the total number of orbitals
    BlockMatrix<Scalar> t;  // the T1-amplitudes as a block matrix, implementing easy operator(i,a) calls


public:
    /*
     *  CONSTRUCTORS
     */

    /**
     *  Construct the T1-amplitudes given their representation as a BlockMatrix.
     * 
     *  @param t                the T1-amplitudes as a block matrix, implementing easy operator(i,a) calls
     *  @param N                the number of occupied orbitals
     *  @param M                the total number of orbitals
     */
    T1Amplitudes(const BlockMatrix<double>& t, const size N, const size_t M) :
        t {t},
        N {N},
        M {M} {}


    /*
     *  OPERATORS
     */

    /**
     *  @param i            an occupied index
     *  @param a            a virtual index
     * 
     *  @return the T1-amplitude corresponding to t_i^a
     */
    double operator()(const size_t i, const size_t a) const { return this->t(i, a); }


    /*
     *  PUBLIC METHODS
     */

    /**
     *  @return the T1-amplitudes as a BlockMatrix
     */
    const BlockMatrix<double>& asBlockMatrix() const { return this->t; }
};


}  // namespace GQCP
