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


#include "Basis/SpinorBasis/Spin.hpp"
#include "ONVBasis/BaseFrozenCoreONVBasis.hpp"
#include "ONVBasis/BaseONVBasis.hpp"
#include "ONVBasis/SpinResolvedONVBasis.hpp"
#include "ONVBasis/SpinUnresolvedFrozenONVBasis.hpp"
#include "Operator/SecondQuantized/USQHamiltonian.hpp"


namespace GQCP {


/**
 *  A frozen spin-resolved ONV basis.
 */
class SpinResolvedFrozenONVBasis: public BaseFrozenCoreONVBasis {
private:
    size_t X;  // number of frozen orbitals/electrons

    SpinUnresolvedFrozenONVBasis frozen_onv_basis_alpha;
    SpinUnresolvedFrozenONVBasis frozen_onv_basis_beta;

    SpinResolvedONVBasis active_onv_basis;  // active (non-frozen) spin-resolved ONV basis containing only the active electrons (N_alpha-X, N_beta-X) and orbitals (K-X)


public:
    // CONSTRUCTORS
    /**
     *  @param K            the total number of orbitals (equal for alpha and beta)
     *  @param N_alpha      the total number of alpha electrons
     *  @param N_beta       the total number of beta electrons
     *  @param X            the number of frozen orbitals and electrons (equal for alpha and beta)
     */
    SpinResolvedFrozenONVBasis(const size_t K, const size_t N_alpha, const size_t N_beta, const size_t X);

    /**
     *  @param onv_basis        (to be frozen) full product spin-resolved ONV basis
     *  @param X                the number of frozen orbitals and electrons (equal for alpha and beta)
     */
    SpinResolvedFrozenONVBasis(const SpinResolvedONVBasis& onv_basis, const size_t X);


    // PUBLIC OVERRIDDEN METHODS

    /**
     *  @return the type of this ONV basis
     */
    ONVBasisType type() const override { return ONVBasisType::SpinResolvedFrozenONVBasis; }


    // PUBLIC METHODS

    /**
     *  @return the active ONV basis
     */
    const SpinResolvedONVBasis& activeONVBasis() const { return this->active_onv_basis; }

    using BaseFrozenCoreONVBasis::evaluateOperatorDense;

    /**
     *  Evaluate the Hamiltonian in a dense matrix
     *
     *  @param usq_hamiltonian                the Hamiltonian expressed in an unrestricted orthonormal basis 
     *  @param diagonal_values                bool to indicate if diagonal values will be calculated
     *
     *  @return the Hamiltonian's evaluation in a dense matrix with the dimensions of this spin-resolved ONV basis
     */
    SquareMatrix<double> evaluateOperatorDense(const USQHamiltonian<double>& usq_hamiltonian, const bool diagonal_values) const;

    using BaseFrozenCoreONVBasis::evaluateOperatorDiagonal;

    /**
     *  Evaluate the diagonal of the Hamiltonian
     *
     *  @param usq_hamiltonian          the Hamiltonian expressed in an unrestricted orthonormal basis 
     *
     *  @return the Hamiltonian's diagonal evaluation in a vector with the dimension of this spin-resolved basis
     */
    VectorX<double> evaluateOperatorDiagonal(const USQHamiltonian<double>& usq_hamiltonian) const;

    /**
     *  Evaluate the unrestricted Hamiltonian in a matrix vector product
     *
     *  @param usq_hamiltonian                the Hamiltonian expressed in an unrestricted orthonormal basis 
     *  @param x                              the vector upon which the evaluation acts 
     *  @param diagonal                       the diagonal evaluated in this spin-resolved basis
     *
     *  @return the Hamiltonian's evaluation in a dense matrix with the dimensions of this spin-resolved basis
     */
    VectorX<double> evaluateOperatorMatrixVectorProduct(const USQHamiltonian<double>& usq_hamiltonian, const VectorX<double>& x, const VectorX<double>& diagonal) const;

    /**
     *  @return the frozen ONV basis for the alpha-spinors
     */
    const SpinUnresolvedFrozenONVBasis& frozenONVBasisAlpha() const { return this->frozen_onv_basis_alpha; }

    /**
     *  @return the frozen ONV basis for the beta-spinors
     */
    const SpinUnresolvedFrozenONVBasis& frozenONVBasisBeta() const { return this->frozen_onv_basis_beta; }

    /**
     *  @return the number of alpha-electrons that this ONV basis describes
     */
    size_t numberOfAlphaElectrons() const { return this->frozen_onv_basis_alpha.numberOfElectrons(); }

    /**
     *  @return the number of beta-electrons that this ONV basis describes
     */
    size_t numberOfBetaElectrons() const { return this->frozen_onv_basis_beta.numberOfElectrons(); }

    /**
     *  @return the number of frozen orbitals that are in this ONV basis
     */
    size_t numberOfFrozenOrbitals() const { return this->X; }
};


}  // namespace GQCP
