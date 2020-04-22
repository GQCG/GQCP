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


/**
 *  Extensions of the <memory> header that haven't been included in C++11
 */


#include <cstddef>
#include <memory>
#include <type_traits>
#include <utility>


/**
 *  A C++11 implementation for std::make_unique. Taken from https://stackoverflow.com/a/17902439.
 */

namespace GQCP {


template <class T>
struct _Unique_if {
    using _Single_object = std::unique_ptr<T>;
};


template <class T>
struct _Unique_if<T[]> {
    using _Unknown_bound = std::unique_ptr<T[]>;
};


template <class T, size_t N>
struct _Unique_if<T[N]> {
    using _Known_bound = void;
};


template <class T, class... Args>
typename _Unique_if<T>::_Single_object make_unique(Args&&... args) {
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}


template <class T>
typename _Unique_if<T>::_Unknown_bound make_unique(size_t n) {
    using U = typename std::remove_extent<T>::type;
    return std::unique_ptr<T>(new U[n]());
}


template <class T, class... Args>
typename _Unique_if<T>::_Known_bound make_unique(Args&&...) = delete;


}  // namespace GQCP
