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


#include "Utilities/Eigen.hpp"


/**
 *  An extension of the Eigen::Array class, with extra operations
 *
 *  @tparam _Scalar     the scalar representation type
 *  @tparam _Rows       the number of rows (int or Dynamic)
 *  @tparam _Cols       the number of columns (int or Dynamic)
 *
 */
namespace GQCP {

template <typename _Scalar = double, int _Rows = Dynamic, int _Cols = Dynamic>
class Array:
    public Eigen::Array<_Scalar, _Rows, _Cols> {


public:
    using Scalar = _Scalar;
    static constexpr auto Rows = _Rows;
    static constexpr auto Cols = _Cols;


    using Self = Array<Scalar, Rows, Cols>;
    using Base = Eigen::Array<Scalar, Rows, Cols>;


public:
    /*
     *  CONSTRUCTORS
     */

    using Eigen::Array<Scalar, Rows, Cols>::Array;  // inhert base constructors


    /**
     *  PUBLIC METHODS
     */

    /**
     *  @return this as a const Eigen::Array
     */
    const Base& Eigen() const { return static_cast<const Base&>(*this); }

    /**
     *  @return this as a non-const Eigen::Matrix
     */
    Base& Eigen() { return static_cast<Base&>(*this); }
};


/*
 *  Convenience typedefs related to Array.
 */

template <typename Scalar>
using ArrayX = Array<Scalar, Dynamic, 1>;

template <typename Scalar>
using ArrayXX = Array<Scalar, Dynamic, Dynamic>;


}  // namespace GQCP
