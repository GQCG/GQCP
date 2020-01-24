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


#include "Mathematical/Representation/SquareMatrix.hpp"

#include <iomanip>


namespace GQCP {


/**
 *  A class that represents a transformation matrix between two orbital bases. The matrix representation of this transformation matrix is such that a new orbital basis b' is found as
 *      b' = b T ,
 *   in which the current orbitals are collected as elements of a row vector b
 * 
 *  @tparam TransformationScalar            the scalar representation of one of the elements of the transformation matrix
 */
template <typename _TransformationScalar>
class TransformationMatrix : public SquareMatrix<_TransformationScalar> {
public:
    using TransformationScalar = _TransformationScalar;


public:
    /*
     *  CONSTRUCTORS
     */

    using SquareMatrix<TransformationScalar>::SquareMatrix;  // inherit base constructors


    /*
     *  PUBLIC METHODS
     */

    /**
     *  In-place 'transform' this transformation matrix such that the resulting transformation matrix describes this and the other transformation together
     */
    void transform(const TransformationMatrix<TransformationScalar>& T) {

        (*this) = (*this) * T;
    }

    /**
     *  @return returns this transformation as to the $VEC group format of GAMESS-US
     */ 
    std::string asGamessUsVecGroup() {
        std::string output = "$VEC";

        std::ostringstream doubleStream;  // Allows conversion of double to the right formatted string
        // Format parameters for $VEC group file
        doubleStream<<std::setprecision(8); 
        doubleStream<<std::scientific;

        for (size_t i = 0; i < this->dimension(); i++) {
            for (size_t j = 0; j < this->dimension(); j++) {

                // Every 5 coefficients is a new line
                if (j % 5 == 0) {
                    output += "\n";
                    output += std::to_string((i+1));
                    output += "  ";
                    output += std::to_string(((j / 5) + 1));
                }

                double x = this->operator()(j ,i);

                // If the double is not negative the space between following entries is padded
                if (x >= 0) {
                    output += " ";
                }

                doubleStream << x;
                std::string double_string = doubleStream.str();

                // GAMESS-US example uses E instead of e
                double_string[double_string.find('e')] = 'E';
                output += double_string;

                // Clear the string
                doubleStream.str("");
            }
        }

        output += "\n$END\n";
        return output;
    }
};


}  // namespace GQCP
