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

        // Extract the current T1- and T2-amplitudes and prepare some variables.
        const auto& sq_hamiltonian = environment.sq_hamiltonian;
        const auto& t1 = environment.t1_amplitudes.back();
        const auto& t2 = environment.t2_amplitudes.back();

        const auto& orbital_space = t1.orbitalSpace();  // assume the orbital spaces are equal for the T1- and T2-amplitudes.
        const auto F = sq_hamiltonian.calculateInactiveFockian(orbital_space).parameters();

        // Update the T1-amplitudes.
        auto t1_updated = t1;
        for (const auto& i : orbital_space.indices(OccupationType::k_occupied)) {
            for (const auto& a : orbital_space.indices(OccupationType::k_virtual)) {

                // Determine the current value for the corresponding T1-amplitude equation, and use it to update the T1-amplitude.
                const auto f_ia = QCModel::CCSD<Scalar>::calculateT1AmplitudeEquation(sq_hamiltonian, t1, t2, i, a);
                t1_updated(i, a) = t1(i, a) + f_ia / (F(i, i) - F(a, a));
            }
        }


        std::cout << "New T1 amplitudes: " << std::endl
                  << t1_updated.asImplicitMatrixSlice().asMatrix() << std::endl
                  << std::endl;

        // Update the T2-amplitudes.
        auto t2_updated = t2;
        for (const auto& i : orbital_space.indices(OccupationType::k_occupied)) {
            for (const auto& j : orbital_space.indices(OccupationType::k_occupied)) {
                for (const auto& a : orbital_space.indices(OccupationType::k_virtual)) {
                    for (const auto& b : orbital_space.indices(OccupationType::k_virtual)) {

                        // Determine the current value for the corresponding T2-amplitude equation, and use it to update the T2-amplitude.
                        const auto f_ijab = QCModel::CCSD<Scalar>::calculateT2AmplitudeEquation(sq_hamiltonian, t1, t2, i, j, a, b);
                        t2_updated(i, j, a, b) = t2(i, j, a, b) + f_ijab / (F(i, i) + F(j, j) - F(a, a) - F(b, b));
                    }
                }
            }
        }


        std::cout << "New T1 amplitudes: " << std::endl
                  << t1_updated.asImplicitMatrixSlice().asMatrix() << std::endl
                  << std::endl;

        // std::cout << "New T2 amplitudes: " << std::endl
        //           << t2_updated.asImplicitRankFourTensorSlice).asTensor() << std::endl
        //           << std::endl;

        // Write the updated amplitudes back to the environment.
        environment.t1_amplitudes.push_back(t1_updated);
        environment.t2_amplitudes.push_back(t2_updated);
    }
};


}  // namespace GQCP
