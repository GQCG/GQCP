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
#include "RDM/TwoRDM.hpp"


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param d    the explicit matrix representation of the 2-RDM
 */
TwoRDM::TwoRDM(const Eigen::Tensor<double, 4>& d) :
    BaseRDM (d.dimensions()[0]),
    d (d)
{
    // Check if the given tensor is 'square'
    auto dims = d.dimensions();
    if ((dims[0] != dims[1]) || (dims[1] != dims[2]) || (dims[2] != dims[3]) ) {
        throw std::invalid_argument("The given tensor should have equal dimensions in every rank.");
    }
}



/*
 *  PUBLIC METHODS
 */

/**
 *  @return the trace of the 2-RDM, i.e. d(p,p,q,q)
 */
double TwoRDM::trace() {
    // TODO: when Eigen3 releases tensor.trace(), use it to implement the reduction

    auto K = static_cast<size_t>(this->d.dimension(1));

    double trace = 0.0;
    for (size_t p = 0; p < K; p++) {
        for (size_t q = 0; q < K; q++) {
            trace += this->d(p,p,q,q);
        }
    }

    return trace;
}


/**
 *  @return a partial contraction of the 2-RDM, where D(p,q) = d(p,q,r,r)
 */
Eigen::MatrixXd TwoRDM::reduce() {
    // TODO: when Eigen3 releases tensor.trace(), use it to implement the reduction

    auto K = static_cast<size_t>(this->d.dimension(1));

    Eigen::MatrixXd D = Eigen::MatrixXd::Zero(K, K);
    for (size_t p = 0; p < K; p++) {
        for (size_t q = 0; q < K; q++) {
            for (size_t r = 0; r < K; r++) {
                D(p,q) += this->d(p,q,r,r);
            }
        }
    }

    return D;
}


}  // namespace GQCP

