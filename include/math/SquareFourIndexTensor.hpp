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
#ifndef SquareFourIndexTensor_hpp
#define SquareFourIndexTensor_hpp


#include "typedefs.hpp"



namespace GQCP {


/**
 *  A square extension of a 4-index tensor class
 *
 *  @tparam Scalar      the scalar type
 */
template<typename Scalar>
class SquareFourIndexTensor: public FourIndexTensor<Scalar> {
public:

    /*
     *  CONSTRUCTORS
     */

    /**
     *  Default constructor
     */
    SquareFourIndexTensor() : FourIndexTensor<Scalar>() {}


    /**
     *  A basic constructor that checks if the given tensor is 'square'
     *
     *  @param tensor       the tensor that should be 'square'
     */
    SquareFourIndexTensor(const FourIndexTensor<Scalar>& tensor) :
        FourIndexTensor<Scalar>(tensor)
    {
        // Check if the given tensor is 'square'
        auto dims = tensor.dimensions();
        if ((dims[0] != dims[1]) || (dims[1] != dims[2]) || (dims[2] != dims[3]) ) {
            throw std::invalid_argument("The given tensor should have equal dimensions in every rank.");
        }
    }


    /*
     *  GETTERS
     */

    size_t get_dim() const {
        return this->dimension(0);  // all tensor dimensions are equal because of the constructor
    }
};


}  // namespace GQCP


#endif /* FourIndexTensor_hpp */
