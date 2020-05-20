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
 *  An iteration step that calculates the new T1 and T2-amplitudes using an update formula from the current T1- and T2-amplitudes.
 * 
 *  @tparam _Scalar             the scalar type that is used the amplitudes
 */
template <typename _Scalar>
class CCSDAmplitudesUpdate:
    public Step<CCSDEnvironment<_Scalar>> {

public:
    using Scalar = _Scalar;
    using Environment = CCSDEnvironment<Scalar>;


public:
    /*
     *  OVERRIDDEN PUBLIC METHODS
     */

    /**
     *  @return a textual description of this algorithmic step
     */
    std::string description() const override {
        return "Calculate the new T1 and T2-amplitudes using an update formula from the current T1- and T2-amplitudes.";
    }

    /**
     *  Calculate the new T1 and T2-amplitudes using an update formula from the current T1- and T2-amplitudes.
     * 
     *  @param environment              the environment that acts as a sort of calculation space
     */
    void execute(Environment& environment) override {

        // Extract the current T1- and T2-amplitudes and intermediates.
        const auto& f = environment.f;
        const auto& V_A = environment.V_A;
        const auto& t1 = environment.t1_amplitudes.back();
        const auto& t2 = environment.t2_amplitudes.back();

        const auto& tau2 = environment.tau2;
        const auto& tau2_tilde = environment.tau2_tilde;

        const auto& F1 = environment.F1;
        const auto& F2 = environment.F2;
        const auto& F3 = environment.F3;

        const auto& W1 = environment.W1;
        const auto& W2 = environment.W2;
        const auto& W3 = environment.W3;

        const auto& orbital_space = t1.orbitalSpace();  // assume the orbital spaces are equal for the T1- and T2-amplitudes.


        // Update the T1-amplitudes.
        auto t1_updated = t1;
        for (const auto& i : orbital_space.indices(OccupationType::k_occupied)) {
            for (const auto& a : orbital_space.indices(OccupationType::k_virtual)) {

                // Determine the current value for the corresponding T1-amplitude equation, and use it to update the T1-amplitude.
                const auto f_ia = QCModel::CCSD<Scalar>::calculateT1AmplitudeEquation(i, a, f, V_A, t1, t2, F1, F2, F3);
                t1_updated(i, a) += f_ia / (f(i, i) - f(a, a));
            }
        }

        // Update the T2-amplitudes.
        auto t2_updated = t2;
        for (const auto& i : orbital_space.indices(OccupationType::k_occupied)) {
            for (const auto& j : orbital_space.indices(OccupationType::k_occupied)) {
                for (const auto& a : orbital_space.indices(OccupationType::k_virtual)) {
                    for (const auto& b : orbital_space.indices(OccupationType::k_virtual)) {

                        // Determine the current value for the corresponding T2-amplitude equation, and use it to update the T2-amplitude.
                        const auto f_ijab = QCModel::CCSD<Scalar>::calculateT2AmplitudeEquation(i, j, a, b, f, V_A, t1, t2, tau2, F1, F2, F3, W1, W2, W3);
                        t2_updated(i, j, a, b) += f_ijab / (f(i, i) + f(j, j) - f(a, a) - f(b, b));
                    }
                }
            }
        }

        // Write the updated amplitudes back to the environment.
        environment.t1_amplitudes.push_back(t1_updated);
        environment.t2_amplitudes.push_back(t2_updated);
    }
};


}  // namespace GQCP
