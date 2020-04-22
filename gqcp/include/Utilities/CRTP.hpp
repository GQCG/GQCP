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


namespace GQCP {


/**
 *  A helper class that helps to implement the CRTP (curiously recurring template pattern) that achieves compile-time polymorphism
 * 
 *  @tparam _Derived            the type of the derived class
 * 
 *  Implementation adapted from (https://www.fluentcpp.com/2017/05/19/crtp-helper/)
 */
template <typename _Derived>
class CRTP {
public:
    using Derived = _Derived;


public:
    /*
     *  PUBLIC METHODS
     */

    /**
     *  @return this as a Derived object (done at compile-time)
     */
    Derived& derived() { return static_cast<Derived&>(*this); }

    /**
     *  @return this as a const Derived object (done at compile-time)
     */
    const Derived& derived() const { return static_cast<const Derived&>(*this); }
};


}  // namespace GQCP
