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


#include "QCModel/HF/Stability/StabilityProperties/StabilityProperties.hpp"


namespace GQCP {


/**
 *  A class that contains the stability information of a GHF wavefunction.
 * 
 */
class GHFStabilityProperties {
public:
    int internal_stability = Stability::unknown;
    int external_stability = Stability::unknown;

public:
    /*
     *  MARK: public methods
     */

    /**
     *  Update the internal stability status of this GHF wavefunction.
     * 
     *  @param status            The new status of the internal stability, this can either be unknown, stable or unstable.
     */
    void updateInternalStability(const Stability status) {
        this->internal_stability = status;
    }


    /**
     *  Update the external stability status of this GHF wavefunction.
     * 
     *  @param status            The new status of the external stability, this can either be unknown, stable or unstable.
     */
    void updateExternalStability(const Stability status) {
        this->external_stability = status;
    }


    /**
     *  Print the stability status of the GHF wavefunction.
     */
    void print() const {

        // First we check the internal stability and print the corresponding message.
        if (this->internal_stability == Stability::unknown) {
            std::cout << "The internal stability of the GHF wavefunction is currently unknown." << std::endl;
        } else if (this->internal_stability == Stability::stable) {
            std::cout << "The GHF wavefunction is internally stable." << std::endl;
        } else if (this->internal_stability == Stability::unstable) {
            std::cout << "The GHF wavefunction is internally unstable." << std::endl;
        }

        // Next we check the external stability and print the corresponding message.
        if (this->external_stability == Stability::unknown) {
            std::cout << "The external stability of the GHF wavefunction is currently unknown." << std::endl;
        } else if (this->external_stability == Stability::stable) {
            std::cout << "The GHF wavefunction is externally stable." << std::endl;
        } else if (this->external_stability == Stability::unstable) {
            std::cout << "The GHF wavefunction is externally unstable." << std::endl;
        }
    }
};

}  // namespace GQCP
