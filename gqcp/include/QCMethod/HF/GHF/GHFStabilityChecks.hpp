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


#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "QCModel/HF/GHF.hpp"
#include "QCModel/HF/Stability/StabilityMatrices/GHFStabilityMatrix.hpp"
#include "Utilities/aliases.hpp"

#include <Eigen/Dense>

#include <algorithm>


namespace GQCP {
namespace QCMethod {


/**
 *  A class that can be used to check the stability conditions of the given GHF wavefunction.
 */
template <typename _Scalar>
class GHFStabilityChecks {
public:
    using Scalar = _Scalar;


public:
    /*
     *  MARK: public methods
     */

    /**
     *  @param ghf_parameters           The GHF model containing the ground state parameters of the wavefunction.
     *  @param gsq_hamiltonian          The second quantized Hamiltonian containing the needed integrals.
     */
    template <typename S = Scalar, typename = IsReal<S>>
    void internalStabilityCheck(GQCP::QCModel::GHF<S>& parameters, const GSQHamiltonian<S>& gsq_hamiltonian) {

        // Get the stability properties of the given GHF wavefunction.
        auto stability = parameters.stabilityProperties();

        // The first step is to calculate the correct stability matrix: This method checks the internal stability of a real valued wavefunction.
        const auto stability_matrix = GHFStabilityMatrix<double>::Internal(parameters, gsq_hamiltonian);

        // Create an eigensolver to diagonalize the stability matrix.
        using MatrixType = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
        Eigen::SelfAdjointEigenSolver<MatrixType> eigensolver {stability_matrix};

        // Calculate the eigenvalues and check whether they are strictly positive or not, we update the stability properties accordingly.
        const auto& eigenvalues = eigensolver.eigenvalues().transpose();
        if (eigenvalues[0] < -1.0e-5) {
            stability.updateInternalStability(Stability::unstable);
        } else {
            stability.updateInternalStability(Stability::stable);
        }

        // Update the stability properties of the GHF wavefunction model.
        parameters.updateStabilityProperties(stability);
    }


    /**
     *  @param ghf_parameters           The GHF model containing the ground state parameters of the wavefunction.
     *  @param gsq_hamiltonian          The second quantized Hamiltonian containing the needed integrals.
     */
    template <typename S = Scalar, typename = IsReal<S>>
    void externalStabilityCheck(GQCP::QCModel::GHF<S>& parameters, const GSQHamiltonian<S>& gsq_hamiltonian) {

        // Get the stability properties of the given GHF wavefunction.
        auto stability = parameters.stabilityProperties();

        // The first step is to calculate the correct stability matrix: This method checks the internal stability of a real valued wavefunction.
        const auto stability_matrix = GHFStabilityMatrix<double>::External(parameters, gsq_hamiltonian);

        // Create an eigensolver to diagonalize the stability matrix.
        using MatrixType = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
        Eigen::SelfAdjointEigenSolver<MatrixType> eigensolver {stability_matrix};

        // Calculate the eigenvalues and check whether they are strictly positive or not, we update the stability properties accordingly.
        const auto& eigenvalues = eigensolver.eigenvalues().transpose();
        if (eigenvalues[0] < -1.0e-5) {
            stability.updateExternalStability(Stability::unstable);
        } else {
            stability.updateExternalStability(Stability::stable);
        }

        // Update the stability properties of the GHF wavefunction model.
        parameters.updateStabilityProperties(stability);
    }
};

}  // namespace QCMethod
}  // namespace GQCP
