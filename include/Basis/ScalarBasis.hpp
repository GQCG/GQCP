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


#include "Basis/ShellSet.hpp"
#include "Mathematical/LinearCombination.hpp"


namespace GQCP {


/**
 *  A class that represents a scalar basis: it provides an interface to obtain basis functions and calculate integrals
 *
 * @tparam _BasisFunction       the type of basis function that this scalar basis 'contains'
 */
template <typename _BasisFunction>
class ScalarBasis {
public:
    using BasisFunction = _BasisFunction;


private:
    ShellSet<BasisFunction> shell_set;  // a collection of shells that can represent this scalar basis


public:

    /*
     *  PUBLIC METHODS
     */

    /**
     *  @return the number of basis functions that 'are' in this scalar basis
     */
    size_t numberOfBasisFunctions() const {
        return this->shell_set.numberOfBasisFunctions();
    }

    /**
     *  @return the basis functions that 'are' in this scalar basis
     */
    std::vector<LinearCombination<double, BasisFunction>> basisFunctions() const {
        static_assert(false, "ScalarBasis::BasisFunctions(): This method has not been implemented yet.");
    }
};


}  // namespace GQCP
