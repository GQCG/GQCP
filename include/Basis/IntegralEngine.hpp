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
#ifndef GQCP_ONEELECTRONINTEGRALENGINE_HPP
#define GQCP_ONEELECTRONINTEGRALENGINE_HPP


#include "Basis/LibintOneElectronIntegralEngine.hpp"
#include "Operator/FirstQuantized/Operator.hpp"
#include "typedefs.hpp"


namespace GQCP {


/**
 *  A class that produces integral engines, much like a factory class
 */
class IntegralEngine {
public:

    /*
     *  LIBINT - ONE-ELECTRON ENGINES
     * 
     *  The following functions create one-electron integral engine using the Libint integral library backend, with the correct template arguments filled in
     */

    // TODO: can be generalized if Operators derive from OneElectronOperator and TwoElectronOperator
    static auto Libint(const OverlapOperator& op) -> LibintOneElectronIntegralEngine<OverlapOperator::Components, product_t<OverlapOperator::Scalar, GTOShell::BasisFunction::Valued>>;

    static auto Libint(const KineticOperator& op) -> LibintOneElectronIntegralEngine<KineticOperator::Components, product_t<KineticOperator::Scalar, GTOShell::BasisFunction::Valued>>;

    static auto Libint(const NuclearAttractionOperator& op) -> LibintOneElectronIntegralEngine<NuclearAttractionOperator::Components, product_t<NuclearAttractionOperator::Scalar, GTOShell::BasisFunction::Valued>>;


    /*
     *  LIBINT - TWO-ELECTRON ENGINES
     * 
     *  The following functions create two-electron integral engine using the Libint integral library backend, with the correct template arguments filled in
     */
    // static auto Libint(const CoulombRepulsionOperator& op) -> LibintTwoElectronIntegralEngine<CoulombRepulsionOperator::Components, product_t<CoulombRepulsionOperator::Scalar, GTOShell::BasisFunction::Valued>>;
};


}  // namespace GQCP


#endif  // GQCP_ONEELECTRONINTEGRALENGINE_HPP
