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


#include "Mathematical/Algorithm/Step.hpp"
#include "QCMethod/CC/CCSDEnvironment.hpp"


namespace GQCP {


/**
 *  An iteration step that calculates the current CCSD electronic correlation energy.
 * 
 *  @tparam _Scalar             the scalar type that is used to represent the amplitudes
 */
template <typename _Scalar>
class CCSDEnergyCalculation:
    public Step<CCSDEnvironment<_Scalar>> {

public:
    using Scalar = _Scalar;
    using Environment = CCSDEnvironment<Scalar>;


public:
    /*
     *  PUBLIC OVERRIDDEN METHODS
     */

    /**
     *  @return a textual description of this algorithmic step
     */
    std::string description() const override {
        return "Calculate the current CCSD electronic correlation energy.";
    }


    /**
     *  Calculate the current CCSD electronic correlation energy.
     * 
     *  @param environment              the environment that acts as a sort of calculation space
     */
    void execute(Environment& environment) override {

        // Prepare some variables.
        const auto& f = environment.f;
        const auto& V_A = environment.V_A;
        const auto& t1 = environment.t1_amplitudes.back();
        const auto& t2 = environment.t2_amplitudes.back();

        // Calculate the current correlation energy and push it to the environment.
        const auto current_correlation_energy = QCModel::CCSD<Scalar>::calculateCorrelationEnergy(f, V_A, t1, t2);
        environment.electronic_energies.push_back(current_correlation_energy);
    }
};


}  // namespace GQCP
