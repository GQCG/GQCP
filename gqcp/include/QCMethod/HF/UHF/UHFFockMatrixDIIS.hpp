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
#include "Mathematical/Optimization/Accelerator/DIIS.hpp"
#include "QCMethod/HF/UHF/UHFFockMatrixDiagonalization.hpp"
#include "QCMethod/HF/UHF/UHFSCFEnvironment.hpp"
#include "QCModel/HF/UHF.hpp"

#include <algorithm>


namespace GQCP {


/**
 *  An iteration step that accelerates the alpha- and beta- Fock matrices (expressed in the scalar/AO basis) based on a DIIS accelerator.
 * 
 *  @tparam _Scalar              the scalar type used to represent the expansion coefficient/elements of the transformation matrix
 */
template <typename _Scalar>
class UHFFockMatrixDIIS:
    public Step<UHFSCFEnvironment<_Scalar>> {

public:
    using Scalar = _Scalar;
    using Environment = UHFSCFEnvironment<Scalar>;


private:
    size_t minimum_subspace_dimension;  // the minimum number of Fock matrices that have to be in the subspace before enabling DIIS
    size_t maximum_subspace_dimension;  // the maximum number of Fock matrices that can be handled by DIIS

    DIIS<Scalar> diis;  // the DIIS accelerator


public:
    /*
     *  CONSTRUCTORS
     */

    /**
     *  @param minimum_subspace_dimension       the minimum number of Fock matrices that have to be in the subspace before enabling DIIS
     *  @param maximum_subspace_dimension       the maximum number of Fock matrices that can be handled by DIIS
     */
    UHFFockMatrixDIIS(const size_t minimum_subspace_dimension = 6, const size_t maximum_subspace_dimension = 6) :
        minimum_subspace_dimension {minimum_subspace_dimension},
        maximum_subspace_dimension {maximum_subspace_dimension} {}


    /*
     *  PUBLIC OVERRIDDEN METHODS
     */

    /**
     *  @return a textual description of this algorithmic step
     */
    std::string description() const override {
        return "Calculate the accelerated alpha- and beta- Fock matrices, and perform a diagonalization step on them.";
    }


    /**
     *  Calculate the accelerated alpha- and beta- Fock matrices, and perform a diagonalization step on them.
     * 
     *  @param environment              the environment that acts as a sort of calculation space
     */
    void execute(Environment& environment) override {

        if (environment.error_vectors_alpha.size() < this->minimum_subspace_dimension) {  // the beta dimension will be the same

            // No acceleration is possible, so calculate the regular Fock matrices and diagonalize them.
            UHFFockMatrixCalculation<Scalar>().execute(environment);
            UHFFockMatrixDiagonalization<Scalar>().execute(environment);
            return;
        }

        // Convert the deques in the environment to vectors that can be accepted by the DIIS accelerator. The total number of elements we can use in DIIS is either the maximum subspace dimension or the number of available error matrices.
        const auto n = std::min(this->maximum_subspace_dimension, environment.error_vectors_alpha.size());
        const std::vector<VectorX<Scalar>> error_vectors_alpha {environment.error_vectors_alpha.end() - n, environment.error_vectors_alpha.end()};  // the n-th last alpha error vectors
        const std::vector<VectorX<Scalar>> error_vectors_beta {environment.error_vectors_beta.end() - n, environment.error_vectors_beta.end()};     // the n-th last beta error vectors

        const std::vector<ScalarUSQOneElectronOperator<Scalar>> fock_matrices {environment.fock_matrices.end() - n, environment.fock_matrices.end()};  // the n-th last Fock matrices

        std::vector<SquareMatrix<Scalar>> alpha_fock_matrices;
        std::vector<SquareMatrix<Scalar>> beta_fock_matrices;

        for (size_t i = 0; i < fock_matrices.size(); i++) {
            alpha_fock_matrices.push_back(fock_matrices[i].alpha().parameters());
            beta_fock_matrices.push_back(fock_matrices[i].beta().parameters());
        }

        // Calculate the accelerated Fock matrices and do a diagonalization step on them.
        const auto F_alpha_accelerated = this->diis.accelerate(alpha_fock_matrices, error_vectors_alpha);
        const auto F_beta_accelerated = this->diis.accelerate(beta_fock_matrices, error_vectors_alpha);

        const ScalarUSQOneElectronOperator<Scalar> F_accelerated {F_alpha_accelerated, F_beta_accelerated};

        environment.fock_matrices.push_back(F_accelerated);  // the diagonalization step can only read from the environment //FIXME

        UHFFockMatrixDiagonalization<Scalar>().execute(environment);

        // The accelerated/extrapolated Fock matrices should not be used in further extrapolation steps, as they are not created from a density matrix.
        environment.fock_matrices.pop_back();
    }
};


}  // namespace GQCP
