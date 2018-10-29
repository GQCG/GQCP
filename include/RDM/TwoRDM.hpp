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
#ifndef GQCP_TWORDM_HPP
#define GQCP_TWORDM_HPP


#include "RDM/BaseRDM.hpp"

#include <unsupported/Eigen/CXX11/Tensor>


namespace GQCP {

/**
 *  A class that holds the tensor representations of a 2RDM
 */
class TwoRDM : public BaseRDM {
private:
    Eigen::Tensor<double, 4> d;


public:
    // CONSTRUCTORS
    explicit TwoRDM(const Eigen::Tensor<double, 4>& d);


    // OPERATORS
    double operator()(size_t p, size_t q, size_t r, size_t s) const { return this->d(p,q,r,s); }


    // GETTERS
    Eigen::Tensor<double, 4> get_matrix_representation() const { return this->d; }


    // PUBLIC METHODS
    /**
     *  @return the trace of the 2-RDM @param: d(p,p,q,q)
     */
    double trace();

    /**
     *  @return a partial contraction of the 2-RDM,
     *  where D(p,q) = d(p,q,r,r)
     */
    Eigen::MatrixXd reduce();
};


}  // namespace GQCP


#endif  // GQCP_TWORDM_HPP
