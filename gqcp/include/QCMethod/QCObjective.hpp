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


#include "Utilities/CRTP.hpp"


namespace GQCP {


/**
 *  A protocol/virtual base class for expressing quantum chemical objectives.
 * 
 *  For a type to conform to the QCObjective protocol, it must implement:
 *      - isSatisfiedWith() that checks if the given model parameters fulfill the objective
 */
template <typename _DerivedQCObjective>
class QCObjective:
    public CRTP<_DerivedQCObjective> {
public:
    using DerivedQCObjective = _DerivedQCObjective;


public:
    // PUBLIC METHODS

    /**
     *  @param model_parameters     the parameters that should be checked against the objective
     * 
     *  @tparam QCModel             the type of the parameters that should be checked against the objective
     * 
     *  @return if the given model parameters fulfill the objective
     */
    template <typename QCModel>
    bool isSatisfiedWith(const QCModel& model_parameters) const {
        return this->derived().isSatisfiedWith(model_parameters);
    }
};


}  // namespace GQCP
