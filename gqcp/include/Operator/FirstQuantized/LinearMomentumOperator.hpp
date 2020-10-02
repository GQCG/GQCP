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


#include "Mathematical/Representation/DenseVectorizer.hpp"
#include "Operator/FirstQuantized/BaseFQOperator.hpp"
#include "Utilities/aliases.hpp"


namespace GQCP {


/**
 *  The (one-electron) vector operator for an electron's linear momentum.
 */
class LinearMomentumOperator:
    public BaseVectorFQOneElectronOperator<complex> {
public:
    /*
     *  MARK: Vectorizer
     */

    // The number of components of the operator.
    static constexpr size_t NumberOfComponents = 3;

    // The 3D vector-vectorizer related to this operator.
    static const VectorVectorizer vectorizer;
};


// Instantiate the static const vectorizer.
const VectorVectorizer LinearMomentumOperator::vectorizer {{LinearMomentumOperator::NumberOfComponents}};


}  // namespace GQCP
