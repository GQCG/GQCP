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


#include "Mathematical/Representation/QCRankFourTensor.hpp"
#include "Processing/RDM/OneRDM.hpp"


namespace GQCP {


/**
 *  A class that represents a 2-RDM
 *
 *  @tparam _Scalar     the scalar type
 */
template <typename _Scalar>
class TwoRDM: public QCRankFourTensor<_Scalar> {
public:
    using Scalar = _Scalar;

    using BaseRepresentation = QCRankFourTensor<Scalar>;


public:
    /*
     *  CONSTRUCTORS
     */

    using QCRankFourTensor<Scalar>::QCRankFourTensor;  // use base constructors


    /*
     *  PUBLIC METHODS
     */

    /**
     *  @return the trace of the 2-RDM, i.e. d(p,p,q,q)
     */
    Scalar trace() const {
        // TODO: when Eigen3 releases tensor.trace(), use it to implement the reduction

        auto K = static_cast<size_t>(this->dimension());

        Scalar trace {};
        for (size_t p = 0; p < K; p++) {
            for (size_t q = 0; q < K; q++) {
                trace += this->operator()(p, p, q, q);
            }
        }

        return trace;
    }


    /**
     *  @return a partial contraction of the 2-RDM, where D(p,q) = d(p,q,r,r)
     */
    OneRDM<Scalar> reduce() const {

        // TODO: when Eigen3 releases tensor.trace(), use it to implement the reduction

        auto K = static_cast<size_t>(this->dimension());

        OneRDM<double> D = OneRDM<double>::Zero(K, K);
        for (size_t p = 0; p < K; p++) {
            for (size_t q = 0; q < K; q++) {
                for (size_t r = 0; r < K; r++) {
                    D(p, q) += this->operator()(p, q, r, r);
                }
            }
        }

        return D;
    }
};


}  // namespace GQCP
