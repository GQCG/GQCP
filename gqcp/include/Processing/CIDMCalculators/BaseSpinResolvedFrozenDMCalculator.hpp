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


#include "DensityMatrix/SpinResolved1DM.hpp"
#include "DensityMatrix/SpinResolved2DM.hpp"


namespace GQCP {


/**
 *  A class capable of calculating 1- and 2-DMs from wave functions expanded in a frozen spin-resolved ONV basis
 */
class BaseSpinResolvedFrozenDMCalculator: public BaseSpinResolvedDMCalculator {
private:
    size_t X;                                                            // number of frozen orbitals/electrons
    std::shared_ptr<BaseSpinResolvedDMCalculator> active_dm_calculator;  // active (non-frozen core) DM builder performing the BaseSpinResolvedDMCalculator interface in the active space with the frozen core CI wave function

public:
    // CONSTRUCTORS

    /**
     *  @param dm_calculator                shared pointer to active (non-frozen core) DM builder
     *  @param X                            the number of frozen orbitals
     */
    BaseSpinResolvedFrozenDMCalculator(const std::shared_ptr<BaseSpinResolvedDMCalculator> dm_calculator, const size_t X);


    // PUBLIC OVERRIDDEN METHODS

    /**
     *  @param x        the coefficient vector representing the wave function
     *
     *  @return all 1-DMs given a coefficient vector
     */
    SpinResolved1DM<double> calculateSpinResolved1DM(const VectorX<double>& x) const override;

    /**
     *  @param x        the coefficient vector representing the wave function
     *
     *  @return all 2-DMs given a coefficient vector
     */
    SpinResolved2DM<double> calculateSpinResolved2DM(const VectorX<double>& x) const override;
};


}  // namespace GQCP
