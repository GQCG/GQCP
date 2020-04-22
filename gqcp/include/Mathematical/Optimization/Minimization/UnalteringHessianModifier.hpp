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


#include "Mathematical/Optimization/Minimization/BaseHessianModifier.hpp"


namespace GQCP {


class UnalteringHessianModifier: public BaseHessianModifier {
public:
    // CONSTRUCTORS
    using BaseHessianModifier::BaseHessianModifier;  // inherit base constructors


    // PUBLIC OVERRIDDEN METHODS

    /**
     *  @param hessian      the current indefinite Hessian
     * 
     *  @return the given Hessian, i.e. do not alter the current hessian
     */
    SquareMatrix<double> operator()(const SquareMatrix<double>& hessian) override;
};


}  // namespace GQCP
