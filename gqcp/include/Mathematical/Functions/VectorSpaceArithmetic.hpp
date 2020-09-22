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


namespace GQCP {


/**
 *  An (abstract) interface for types that support vector space arithmetic.
 */
template <typename T, typename Scalar>
class VectorSpaceArithmetic {
public:
    /*
     *  MARK: Pure virtual functions
     */

    /**
     *  Addition-assignment, to be implemented in derived classes.
     */
    virtual T& operator+=(const T& rhs) = 0;

    /**
     *  Scalar multiplication-assignment, to be implemented in derived classes.
     */
    virtual T& operator*=(const Scalar& a) = 0;


    /*
     *  MARK: Canonical implementations of operators
     */

    /**
     *  Addition, canonically implemented using addition-assignment.
     */
    friend T operator+(T lhs, const T& rhs) {
        lhs += rhs;
        return lhs;
    }


    /**
     *  Subtraction-assignment, canonically implemented using addition-assignment with the negation.
     */
    T& operator-=(const T& rhs) {
        *this += (-rhs);
        return static_cast<T&>(*this);
    }


    /**
     *  Subtraction, canonically implemented using subtraction-assignment.
     */
    friend T operator-(T lhs, const T& rhs) {
        lhs -= rhs;
        return lhs;
    }


    /**
     *  Scalar multiplication, canonically implemented using scalar multiplication-assignment.
     */
    friend T operator*(const Scalar& a, T rhs) {
        return rhs *= a;
    }

    /**
     *  The commutative version of the previous scalar multiplication.
     */
    friend T operator*(T lhs, const Scalar& a) {
        return lhs *= a;
    }


    /**
     *  Negation, canonically implemented as scalar multiplication by (-1.0).
     */
    T operator-() const {
        return Scalar {-1.0} * static_cast<const T&>(*this);
    };
};


}  // namespace GQCP
