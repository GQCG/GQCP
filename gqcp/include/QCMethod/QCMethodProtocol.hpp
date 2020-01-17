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


#include "QCMethod/QCObjective.hpp"
#include "QCMethod/QCStructure.hpp"


namespace GQCP {


/**
 *  A protocol/virtual base class for quantum chemical methods
 * 
 *  @tparam _QCModel                the type of the quantum chemical model this quantum chemical method is trying to solve
 *  @tparam _DerivedQCMethod        the type of the class that derives from this base class
 */
template <typename _QCModel, typename _DerivedQCMethod>
class QCMethodProtocol: public CRTP<_DerivedQCMethod> {
public:
    using QCModel = _QCModel;
    using DerivedQCMethod = _DerivedQCMethod;


public:

    // PUBLIC METHODS

    /**
     *  Optimize the electronic structure model: find the parameters that are the solutions to the quantum chemical method's objective
     * 
     *  @tparam QCObjective         the type of the objective
     *  @tparam Solver              the type of the solver
     * 
     *  @param objective            the objective that should be fulfilled in order to consider the model's parameters as 'optimal'
     *  @param solver               the solver that will try to optimize the parameters
     */
    template <typename QCObjective, typename Solver>
    QCStructure<QCModel> optimize(const QCObjective& objective, Solver& solver);
};



}  // namespace GQCP
