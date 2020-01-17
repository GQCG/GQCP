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


#include "QCMethod/QCMethodProtocol.hpp"
#include "QCMethod/QCObjective.hpp"
#include "QCModel/HF/RHF.hpp"


namespace GQCP {
namespace QCMethod {


// /**
//  *  The restricted Hartree-Fock quantum chemical method
//  * 
//  *  @tparam _ExpansionScalar_       the type of scalar that is used for the expansion of the spatial orbitals in their underlying scalar basis
//  */
// template <typename _ExpansionScalar>
// class RHF: public GQCP::QCMethodProtocol<QCModel::RHF<_ExpansionScalar>, QCMethod::RHF<_ExpansionScalar>> {
// public:
//     using ExpansionScalar = _ExpansionScalar;

// public:

//     /*
//      *  PUBLIC METHODS
//      */

//     /**
//      *  Optimize the electronic structure model: find the parameters that are the solutions to the quantum chemical method's objective
//      * 
//      *  @tparam QCObjective         the type of the objective
//      *  @tparam Solver              the type of the solver
//      * 
//      *  @param objective            the objective that should be fulfilled in order to consider the model's parameters as 'optimal'
//      *  @param solver               the solver that will try to optimize the parameters
//      */
//     template <typename QCObjective, typename Solver>
//     QCStructure<QCModel::RHF<ExpansionScalar>> optimize(const QCObjective& objective, Solver& solver) {
//         solver.solve();
//         return QCStructure<QCModel::RHF<ExpansionScalar>>(solver.solution());
//     }
// };


}  // namespace QCMethod
}  // namespace GQCP
