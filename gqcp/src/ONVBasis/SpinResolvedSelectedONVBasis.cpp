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

#include "ONVBasis/SpinResolvedSelectedONVBasis.hpp"

#include <boost/dynamic_bitset.hpp>
#include <boost/math/special_functions.hpp>
#include <boost/numeric/conversion/converter.hpp>


namespace GQCP {


/*
 *  MARK: Constructors
 */

/**
 *  Construct an empty spin-resolved selected ONV basis.
 *
 *  @param K            The number of spin-orbitals (equal for alpha and beta).
 *  @param N_alpha      The number of alpha electrons, i.e. the number of occupied alpha spin-orbitals.
 *  @param N_beta       The number of beta electrons, i.e. the number of occupied beta spin-orbitals.
 */
SpinResolvedSelectedONVBasis::SpinResolvedSelectedONVBasis(const size_t K, const size_t N_alpha, const size_t N_beta) :
    K {K},
    N_alpha {N_alpha},
    N_beta {N_beta} {}


/**
 *  Generate a `SpinResolvedSelectedONVBasis` from a seniority-zero ONV basis.
 *
 *  @param onv_basis        The seniority-zero ONV basis.
 */
SpinResolvedSelectedONVBasis::SpinResolvedSelectedONVBasis(const SeniorityZeroONVBasis& onv_basis) :
    SpinResolvedSelectedONVBasis(onv_basis.numberOfSpatialOrbitals(), onv_basis.numberOfElectronPairs(), onv_basis.numberOfElectronPairs()) {

    // Prepare some variables.
    const auto dimension = onv_basis.dimension();
    const auto proxy_onv_basis = onv_basis.proxy();

    // Iterate over the seniority-zero ONV basis and add all ONVs as doubly-occupied ONVs.
    std::vector<SpinResolvedONV> onvs;
    SpinUnresolvedONV onv = proxy_onv_basis.constructONVFromAddress(0);
    for (size_t I = 0; I < dimension; I++) {  // I iterates over all addresses of the doubly-occupied ONVs

        onvs.emplace_back(onv, onv);

        if (I < dimension - 1) {  // prevent the last permutation from occurring
            proxy_onv_basis.transformONVToNextPermutation(onv);
        }
    }

    this->onvs = onvs;
}


/**
 *  Generate a `SpinResolvedSelectedONVBasis` from a full spin-resolved ONV basis.
 *
 *  @param onv_basis        The full spin-resolved ONV basis.
 */
SpinResolvedSelectedONVBasis::SpinResolvedSelectedONVBasis(const SpinResolvedONVBasis& onv_basis) :
    SpinResolvedSelectedONVBasis(onv_basis.alpha().numberOfOrbitals(), onv_basis.alpha().numberOfElectrons(), onv_basis.beta().numberOfElectrons()) {

    std::vector<SpinResolvedONV> onvs;

    const SpinUnresolvedONVBasis& onv_basis_alpha = onv_basis.alpha();
    const SpinUnresolvedONVBasis& onv_basis_beta = onv_basis.beta();

    auto dim_alpha = onv_basis_alpha.dimension();
    auto dim_beta = onv_basis_beta.dimension();

    SpinUnresolvedONV alpha = onv_basis_alpha.constructONVFromAddress(0);
    for (size_t I_alpha = 0; I_alpha < dim_alpha; I_alpha++) {

        SpinUnresolvedONV beta = onv_basis_beta.constructONVFromAddress(0);
        for (size_t I_beta = 0; I_beta < dim_beta; I_beta++) {

            onvs.push_back(SpinResolvedONV {alpha, beta});

            if (I_beta < dim_beta - 1) {  // prevent the last permutation from occurring
                onv_basis_beta.transformONVToNextPermutation(beta);
            }
        }
        if (I_alpha < dim_alpha - 1) {  // prevent the last permutation from occurring
            onv_basis_alpha.transformONVToNextPermutation(alpha);
        }
    }
    this->onvs = onvs;
}


/*
 *  MARK: Named constructors
 */

/**
 *  Create a `SpinResolvedSelectedONVBasis` for a CI singles calculation, using the HF determinant as a reference.
 * 
 *  @param K            The number of spin-orbitals (equal for alpha and beta).
 *  @param N_alpha      The number of alpha electrons, i.e. the number of occupied alpha spin-orbitals.
 *  @param N_beta       The number of beta electrons, i.e. the number of occupied beta spin-orbitals.
 * 
 *  @returen A CI singles-equivalent `SpinResolvedSelectedONVBasis`.
 */
SpinResolvedSelectedONVBasis SpinResolvedSelectedONVBasis::CIS(const size_t K, const size_t N_alpha, const size_t N_beta) {

    const auto V_alpha = K - N_alpha;          // The number of alpha virtual orbitals.
    const auto dim_alpha = N_alpha * V_alpha;  // The number of alpha excitations.

    const auto V_beta = K - N_beta;         // The number of beta virtual orbitals.
    const auto dim_beta = N_beta * V_beta;  // The number of beta excitations.


    SpinResolvedSelectedONVBasis onv_basis {K, N_alpha, N_beta};

    std::vector<SpinResolvedONV> onvs;
    onvs.reserve(dim_alpha * dim_beta);

    auto reference = SpinResolvedONV::UHF(K, N_alpha, N_beta);
    auto alpha_reference = reference.onv(Spin::alpha);
    auto beta_reference = reference.onv(Spin::beta);

    const auto alpha_orbital_space = alpha_reference.orbitalSpace();
    const auto beta_orbital_space = beta_reference.orbitalSpace();

    onvs.push_back(reference);
    // Generate the alpha-alpha-excitations.
    std::cout << "AA" << std::endl;
    for (const auto& i_alpha : alpha_orbital_space.indices(OccupationType::k_occupied)) {
        std::cout << "i: " << i_alpha << std::endl;
        for (const auto& a_alpha : alpha_orbital_space.indices(OccupationType::k_virtual)) {
            std::cout << "a: " << a_alpha << std::endl;
            auto alpha_part = alpha_reference;
            auto beta_part = beta_reference;

            alpha_part.annihilate(i_alpha);
            alpha_part.create(a_alpha);

            onvs.emplace_back(alpha_part, beta_part);
        }
    }

    // Generate the beta-beta-excitations.
    std::cout << "BB" << std::endl;
    for (const auto& i_beta : beta_orbital_space.indices(OccupationType::k_occupied)) {
        std::cout << "i: " << i_beta << std::endl;
        for (const auto& a_beta : beta_orbital_space.indices(OccupationType::k_virtual)) {
            std::cout << "a: " << a_beta << std::endl;
            auto alpha_part = alpha_reference;
            auto beta_part = beta_reference;

            beta_part.annihilate(i_beta);
            beta_part.create(a_beta);

            onvs.emplace_back(alpha_part, beta_part);
        }
    }


    // Generate the alpha-beta excitations.
    std::cout << "AB" << std::endl;
    for (const auto& i_alpha : alpha_orbital_space.indices(OccupationType::k_occupied)) {
        std::cout << "i: " << i_alpha << std::endl;
        for (const auto& a_beta : beta_orbital_space.indices(OccupationType::k_virtual)) {
            std::cout << "a: " << a_beta << std::endl;
            auto alpha_part = alpha_reference;
            auto beta_part = beta_reference;

            alpha_part.annihilate(i_alpha);
            beta_part.create(a_beta);

            onvs.emplace_back(alpha_part, beta_part);
        }
    }


    // Generate the beta-alpha excitations.
    std::cout << "BA" << std::endl;
    for (const auto& i_beta : beta_orbital_space.indices(OccupationType::k_occupied)) {
        std::cout << "i: " << i_beta << std::endl;
        for (const auto& a_alpha : alpha_orbital_space.indices(OccupationType::k_virtual)) {
            std::cout << "a: " << a_alpha << std::endl;
            auto alpha_part = alpha_reference;
            auto beta_part = beta_reference;

            alpha_part.create(a_alpha);
            beta_part.annihilate(i_beta);

            onvs.emplace_back(alpha_part, beta_part);
        }
    }


    onv_basis.expandWith(onvs);
    std::cout << "number of ONVS: " << onv_basis.dimension() << std::endl;

    for (size_t i = 0; i < onv_basis.dimension(); i++) {
        std::cout << onv_basis.onvWithIndex(i).asString() << std::endl;
    }


    return onv_basis;
}


/*
 *  MARK: Modifying
 */

/**
 *  Expand this ONV basis with the given spin-resolved ONV.
 * 
 *  @param onv          The ONV that should be included in this ONV basis.
 */
void SpinResolvedSelectedONVBasis::expandWith(const SpinResolvedONV& onv) {

    if ((onv.onv(Spin::alpha).numberOfElectrons() != this->numberOfAlphaElectrons()) || (onv.onv(Spin::beta).numberOfElectrons() != this->numberOfBetaElectrons())) {
        throw std::invalid_argument("SpinResolvedSelectedONVBasis::expandWith(const SpinResolvedONV&): The given ONV's number of electrons is not compatible with the number of electrons for this ONV basis.");
    }

    if ((onv.onv(Spin::alpha).numberOfSpinors() != this->numberOfOrbitals()) || (onv.onv(Spin::beta).numberOfSpinors() != this->numberOfOrbitals())) {
        throw std::invalid_argument("SpinResolvedSelectedONVBasis::expandWith(const SpinResolvedONV&): The given ONV's number of orbitals is not compatible with the number of orbitals for this ONV basis.");
    }

    this->onvs.push_back(onv);
}


/**
 *  Expand this ONV basis with the given spin-resolved ONVs.
 * 
 *  @param onvs         The ONVs that should be included in this ONV basis.
 */
void SpinResolvedSelectedONVBasis::expandWith(const std::vector<SpinResolvedONV>& onvs) {

    for (const auto& onv : onvs) {
        this->expandWith(onv);
    }
}


/*
 *  MARK: Dense restricted operator evaluations
 */

/**
 *  Calculate the dense matrix representation of a restricted one-electron operator in this ONV basis.
 *
 *  @param f                A restricted one-electron operator expressed in an orthonormal orbital basis.
 *
 *  @return A dense matrix represention of the one-electron operator.
 */
SquareMatrix<double> SpinResolvedSelectedONVBasis::evaluateOperatorDense(const ScalarRSQOneElectronOperator<double>& f) const {

    // By delegating the actual implementation of this method to its unrestricted counterpart, we avoid code duplication for the restricted part.
    // This does not affect performance significantly, because the bottleneck will always be the double iteration over the whole ONV basis.
    const auto f_unrestricted = ScalarUSQOneElectronOperator<double>::FromRestricted(f);
    return this->evaluateOperatorDense(f_unrestricted);
}


/**
 *  Calculate the dense matrix representation of a restricted two-electron operator in this ONV basis.
 *
 *  @param g                A restricted two-electron operator expressed in an orthonormal orbital basis.
 *
 *  @return A dense matrix represention of the two-electron operator.
 */
SquareMatrix<double> SpinResolvedSelectedONVBasis::evaluateOperatorDense(const ScalarRSQTwoElectronOperator<double>& g) const {

    // By delegating the actual implementation of this method to its unrestricted counterpart, we avoid code duplication for the restricted part.
    // This does not affect performance significantly, because the bottleneck will always be the double iteration over the whole ONV basis.
    // Furthermore, we can use the `USQHamiltonian`'s general evaluation function, because even adding zero-valued one-electron operators won't have an impact. This would be different if we would split up the evaluation in one- and two-electron operator evaluations, which would require two times the double iterations over the whole ONV basis.
    const auto zero = ScalarUSQOneElectronOperator<double>::Zero(g.numberOfOrbitals());
    const auto g_unrestricted = ScalarUSQTwoElectronOperator<double>::FromRestricted(g);
    const USQHamiltonian<double> hamiltonian {zero, g_unrestricted};

    return this->evaluateOperatorDense(hamiltonian);
}


/**
 *  Calculate the dense matrix representation of a restricted Hamiltonian in this ONV basis.
 *
 *  @param hamiltonian      A restricted Hamiltonian expressed in an orthonormal orbital basis.
 *
 *  @return A dense matrix represention of the Hamiltonian.
 */
SquareMatrix<double> SpinResolvedSelectedONVBasis::evaluateOperatorDense(const RSQHamiltonian<double>& hamiltonian) const {

    // By delegating the actual implementation of this method to its unrestricted counterpart, we avoid code duplication for the restricted part.
    // This does not affect performance significantly, because the bottleneck will always be the double iteration over the whole ONV basis.
    const auto h_unrestricted = ScalarUSQOneElectronOperator<double>::FromRestricted(hamiltonian.core());
    const auto g_unrestricted = ScalarUSQTwoElectronOperator<double>::FromRestricted(hamiltonian.twoElectron());
    const USQHamiltonian<double> unrestricted_hamiltonian {h_unrestricted, g_unrestricted};

    return this->evaluateOperatorDense(unrestricted_hamiltonian);
}


/*
 *  MARK: Sparse restricted operator evaluations
 */

/**
 *  Calculate the sparse matrix representation of a restricted one-electron operator in this ONV basis.
 *
 *  @param f                A restricted one-electron operator expressed in an orthonormal orbital basis.
 *
 *  @return A sparse matrix represention of the one-electron operator.
 */
Eigen::SparseMatrix<double> SpinResolvedSelectedONVBasis::evaluateOperatorSparse(const ScalarRSQOneElectronOperator<double>& f) const {

    if (f.numberOfOrbitals() != this->numberOfOrbitals()) {
        throw std::invalid_argument("SpinResolvedSelectedONVBasis::evaluateOperatorSparse(const ScalarRSQOneElectronOperator<double>&): The number of orbitals of the ONV basis and the operator are incompatible.");
    }

    // Initialize a container for the sparse matrix representation, and reserve an appropriate amount of memory for it.
    MatrixRepresentationEvaluationContainer<Eigen::SparseMatrix<double>> container {this->dimension()};
    size_t memory = this->dimension() + this->dimension() * this->K * (this->N_alpha + this->N_beta);
    container.reserve(memory);

    // Evaluate the one-electron operator (as an unrestricted operator) and add the evaluations to the sparse matrix representation.
    const auto f_unrestricted = ScalarUSQOneElectronOperator<double>::FromRestricted(f);
    this->evaluate<Eigen::SparseMatrix<double>>(f_unrestricted, container);

    // Finalize the creation of the sparse matrix and return the result.
    container.addToMatrix();
    return container.evaluation();
}


/**
 *  Calculate the sparse matrix representation of a restricted two-electron operator in this ONV basis.
 *
 *  @param g                A restricted two-electron operator expressed in an orthonormal orbital basis.
 *
 *  @return A sparse matrix represention of the two-electron operator.
 */
Eigen::SparseMatrix<double> SpinResolvedSelectedONVBasis::evaluateOperatorSparse(const ScalarRSQTwoElectronOperator<double>& g) const {

    if (g.numberOfOrbitals() != this->numberOfOrbitals()) {
        throw std::invalid_argument("SpinResolvedSelectedONVBasis::evaluateOperatorSparse(const ScalarRSQTwoElectronOperator<double>&): The number of orbitals of the ONV basis and the operator are incompatible.");
    }

    // Initialize a container for the sparse matrix representation, and reserve an appropriate amount of memory for it.
    MatrixRepresentationEvaluationContainer<Eigen::SparseMatrix<double>> container {this->dimension()};

    size_t memory = this->dimension() + this->dimension() * this->K * this->K * (this->N_alpha + this->N_beta) * (this->N_alpha + this->N_beta);
    container.reserve(memory);

    // Use the `USQHamiltonian`'s general evaluation function, because even adding zero-valued one-electron operators won't have an impact. This would be different if we would split up the evaluation in one- and two-electron operator evaluations, which would require two times the double iterations over the whole ONV basis.
    const auto zero = ScalarUSQOneElectronOperator<double>::Zero(g.numberOfOrbitals());
    const auto g_unrestricted = ScalarUSQTwoElectronOperator<double>::FromRestricted(g);
    const USQHamiltonian<double> hamiltonian {zero, g_unrestricted};

    this->evaluate<Eigen::SparseMatrix<double>>(hamiltonian, container);

    // Finalize the creation of the sparse matrix and return the result.
    container.addToMatrix();
    return container.evaluation();
}


/**
 *  Calculate the sparse matrix representation of a restricted Hamiltonian in this ONV basis.
 *
 *  @param hamiltonian      A restricted Hamiltonian expressed in an orthonormal orbital basis.
 *
 *  @return A sparse matrix represention of the Hamiltonian.
 */
Eigen::SparseMatrix<double> SpinResolvedSelectedONVBasis::evaluateOperatorSparse(const RSQHamiltonian<double>& hamiltonian) const {

    // Delegate the implementation to the unrestricted evaluation.
    const auto h_unrestricted = ScalarUSQOneElectronOperator<double>::FromRestricted(hamiltonian.core());
    const auto g_unrestricted = ScalarUSQTwoElectronOperator<double>::FromRestricted(hamiltonian.twoElectron());
    const USQHamiltonian<double> unrestricted_hamiltonian {h_unrestricted, g_unrestricted};

    return this->evaluateOperatorSparse(unrestricted_hamiltonian);
}


/*
 *  MARK: Diagonal restricted operator evaluations
 */

/**
 *  Calculate the diagonal of the matrix representation of a restricted one-electron operator in this ONV basis.
 *
 *  @param f_op             A restricted one-electron operator expressed in an orthonormal orbital basis.
 *
 *  @return The diagonal of the dense matrix represention of the one-electron operator.
 */
VectorX<double> SpinResolvedSelectedONVBasis::evaluateOperatorDiagonal(const ScalarRSQOneElectronOperator<double>& f_op) const {

    const auto K = f_op.numberOfOrbitals();
    if (K != this->numberOfOrbitals()) {
        throw std::invalid_argument("SpinResolvedSelectedONVBasis::evaluateOperatorDiagonal(const ScalarRSQOneElectronOperator<double>&): The number of orbitals of this ONV basis and the given operator are incompatible.");
    }

    // Prepare some variables.
    const auto dim = this->dimension();
    const auto& f = f_op.parameters();
    VectorX<double> diagonal = VectorX<double>::Zero(dim);

    for (size_t I = 0; I < dim; I++) {  // I loops over the addresses of alpha onvs
        SpinResolvedONV onv_I = this->onvWithIndex(I);
        SpinUnresolvedONV alpha_I = onv_I.onv(Spin::alpha);
        SpinUnresolvedONV beta_I = onv_I.onv(Spin::beta);

        for (size_t p = 0; p < K; p++) {
            if (alpha_I.isOccupied(p)) {
                diagonal(I) += f(p, p);
            }

            if (beta_I.isOccupied(p)) {
                diagonal(I) += f(p, p);
            }
        }

    }  // I loop

    return diagonal;
}


/**
 *  Calculate the diagonal of the matrix representation of a restricted two-electron operator in this ONV basis.
 *
 *  @param g_op             A restricted two-electron operator expressed in an orthonormal orbital basis.
 *
 *  @return The diagonal of the dense matrix represention of the two-electron operator.
 */
VectorX<double> SpinResolvedSelectedONVBasis::evaluateOperatorDiagonal(const ScalarRSQTwoElectronOperator<double>& g_op) const {

    const auto K = g_op.numberOfOrbitals();
    if (K != this->numberOfOrbitals()) {
        throw std::invalid_argument("SpinResolvedSelectedONVBasis::evaluateOperatorDiagonal(const ScalarRSQTwoElectronOperator<double>&): The number of orbitals of this ONV basis and the given operator are incompatible.");
    }

    // Prepare some variables.
    const auto dim = this->dimension();
    const auto& g = g_op.parameters();
    VectorX<double> diagonal = VectorX<double>::Zero(dim);

    for (size_t I = 0; I < dim; I++) {  // I loops over addresses of all ONVs
        SpinResolvedONV onv_I = this->onvWithIndex(I);
        SpinUnresolvedONV alpha_I = onv_I.onv(Spin::alpha);
        SpinUnresolvedONV beta_I = onv_I.onv(Spin::beta);

        for (size_t p = 0; p < K; p++) {
            if (alpha_I.isOccupied(p)) {
                for (size_t q = 0; q < K; q++) {

                    if (p != q) {  // can't create/annihilate the same orbital twice
                        if (alpha_I.isOccupied(q)) {
                            diagonal(I) += 0.5 * g(p, p, q, q);
                            diagonal(I) -= 0.5 * g(p, q, q, p);
                        }
                    }

                    if (beta_I.isOccupied(q)) {
                        diagonal(I) += 0.5 * g(p, p, q, q);
                    }
                }  // loop over q
            }

            if (beta_I.isOccupied(p)) {
                for (size_t q = 0; q < K; q++) {

                    if (p != q) {  // can't create/annihilate the same orbital twice
                        if (beta_I.isOccupied(q)) {
                            diagonal(I) += 0.5 * g(p, p, q, q);
                            diagonal(I) -= 0.5 * g(p, q, q, p);
                        }
                    }

                    if (alpha_I.isOccupied(q)) {
                        diagonal(I) += 0.5 * g(p, p, q, q);
                    }
                }  // loop over q
            }
        }  // loop over q

    }  // I loop

    return diagonal;
};

/**
 *  Calculate the diagonal of the dense matrix representation of a restricted Hamiltonian in this ONV basis.
 *
 *  @param hamiltonian      A restricted Hamiltonian expressed in an orthonormal orbital basis.
 *
 *  @return The diagonal of the dense matrix represention of the Hamiltonian.
 */
VectorX<double> SpinResolvedSelectedONVBasis::evaluateOperatorDiagonal(const RSQHamiltonian<double>& hamiltonian) const {

    return this->evaluateOperatorDiagonal(hamiltonian.core()) + this->evaluateOperatorDiagonal(hamiltonian.twoElectron());
}


/*
 *  MARK: Dense unrestricted operator evaluations
 */

/**
 *  Calculate the dense matrix representation of an unrestricted one-electron operator in this ONV basis.
 *
 *  @param f                An unrestricted one-electron operator expressed in an orthonormal orbital basis.
 *
 *  @return A dense matrix represention of the one-electron operator.
 */
SquareMatrix<double> SpinResolvedSelectedONVBasis::evaluateOperatorDense(const ScalarUSQOneElectronOperator<double>& f) const {

    const auto K = this->numberOfOrbitals();
    if ((f.alpha().numberOfOrbitals() != K) || (f.beta().numberOfOrbitals() != K)) {
        throw std::invalid_argument("SpinResolvedSelectedONVBasis::evaluateOperatorDense(const ScalarUSQOneElectronOperator<double>&): The number of orbitals of this ONV basis and the given one-electron operator are incompatible.");
    }

    // Initialize a container for the dense matrix representation, and fill it with the general evaluation function.
    MatrixRepresentationEvaluationContainer<SquareMatrix<double>> container {this->dimension()};
    this->evaluate<SquareMatrix<double>>(f, container);

    return container.evaluation();
}


/**
 *  Calculate the dense matrix representation of an unrestricted Hamiltonian in this ONV basis.
 *
 *  @param hamiltonian      An unrestricted Hamiltonian expressed in an orthonormal orbital basis.
 *
 *  @return A dense matrix represention of the Hamiltonian.
 */
SquareMatrix<double> SpinResolvedSelectedONVBasis::evaluateOperatorDense(const USQHamiltonian<double>& hamiltonian) const {

    if (hamiltonian.numberOfOrbitals() != this->numberOfOrbitals()) {
        throw std::invalid_argument("SpinResolvedSelectedONVBasis::evaluateOperatorDense(const USQHamiltonian<double>&): The number of orbitals of this ONV basis and the given Hamiltonian are incompatible.");
    }

    // Initialize a container for the dense matrix representation, and fill it with the general evaluation function.
    MatrixRepresentationEvaluationContainer<SquareMatrix<double>> container {this->dimension()};
    this->evaluate<SquareMatrix<double>>(hamiltonian, container);

    return container.evaluation();
}


/*
 *  MARK: Diagonal unrestricted operator evaluations
 */

/**
 *  Calculate the diagonal of the dense matrix representation of an unrestricted Hamiltonian in this ONV basis.
 *
 *  @param hamiltonian      An unrestricted Hamiltonian expressed in an orthonormal orbital basis.
 *
 *  @return The diagonal of the dense matrix represention of the Hamiltonian.
 */
VectorX<double> SpinResolvedSelectedONVBasis::evaluateOperatorDiagonal(const USQHamiltonian<double>& hamiltonian) const {

    const auto K = hamiltonian.numberOfOrbitals();
    if (K != this->numberOfOrbitals()) {
        throw std::invalid_argument("SpinResolvedSelectedONVBasis::evaluateOperatorDiagonal(USQHamiltonian<double>): Basis functions of this ONV basis and the operator are incompatible.");
    }

    // Prepare some variables.
    const auto dim = this->dimension();
    const auto& h_a = hamiltonian.core().alpha().parameters();
    const auto& g_aa = hamiltonian.twoElectron().alphaAlpha().parameters();
    const auto& h_b = hamiltonian.core().beta().parameters();
    const auto& g_bb = hamiltonian.twoElectron().betaBeta().parameters();

    // For the mixed two-electron integrals g_ab and g_ba, we can use the following relation: g_ab(pqrs) = g_ba(rspq) and proceed to only work with g_ab.
    const auto& g_ab = hamiltonian.twoElectron().alphaBeta().parameters();


    VectorX<double> diagonal = VectorX<double>::Zero(dim);
    for (size_t I = 0; I < dim; I++) {  // Ia loops over addresses of alpha onvs
        SpinResolvedONV onv_I = this->onvWithIndex(I);
        SpinUnresolvedONV alpha_I = onv_I.onv(Spin::alpha);
        SpinUnresolvedONV beta_I = onv_I.onv(Spin::beta);

        for (size_t p = 0; p < K; p++) {
            if (alpha_I.isOccupied(p)) {

                diagonal(I) += h_a(p, p);

                for (size_t q = 0; q < K; q++) {

                    if (p != q) {  // can't create/annihilate the same orbital twice
                        if (alpha_I.isOccupied(q)) {
                            diagonal(I) += 0.5 * g_aa(p, p, q, q);
                            diagonal(I) -= 0.5 * g_aa(p, q, q, p);
                        }
                    }

                    if (beta_I.isOccupied(q)) {
                        diagonal(I) += 0.5 * g_ab(p, p, q, q);
                    }
                }  // loop over q
            }

            if (beta_I.isOccupied(p)) {

                diagonal(I) += h_b(p, p);

                for (size_t q = 0; q < K; q++) {

                    if (p != q) {  // can't create/annihilate the same orbital twice
                        if (beta_I.isOccupied(q)) {
                            diagonal(I) += 0.5 * g_bb(p, p, q, q);
                            diagonal(I) -= 0.5 * g_bb(p, q, q, p);
                        }
                    }

                    if (alpha_I.isOccupied(q)) {
                        diagonal(I) += 0.5 * g_ab(q, q, p, p);
                    }
                }  // loop over q
            }
        }  // loop over q

    }  // alpha address (Ia) loop

    return diagonal;
}


/*
 *  MARK: Restricted matrix-vector product evaluations
 */

/**
 *  Calculate the matrix-vector product of (the matrix representation of) a one-electron operator with the given coefficient vector.
 *
 *  @param f                A restricted one-electron operator expressed in an orthonormal orbital basis.
 *  @param x                The coefficient vector of a linear expansion.
 *
 *  @return The coefficient vector of the linear expansion after being acted on with the given (matrix representation of) the one-electron operator.
 */
VectorX<double> SpinResolvedSelectedONVBasis::evaluateOperatorMatrixVectorProduct(const ScalarRSQOneElectronOperator<double>& f, const VectorX<double>& x) const {

    if (f.numberOfOrbitals() != this->numberOfOrbitals()) {
        throw std::invalid_argument("SpinResolvedSelectedONVBasis::evaluateOperatorMatrixVectorProduct(const ScalarRSQOneElectronOperator<double>&, const VectorX<double>&): The number of orbitals of this ONV basis and the operator are incompatible.");
    }

    // Convert the restricted operator to an unrestricted operator, and use the general unrestricted one-electron operator evaluation.
    const auto f_unrestricted = ScalarUSQOneElectronOperator<double>::FromRestricted(f);

    // Initialize a container for the matrix-vector product, and fill it with the general evaluation function.
    MatrixRepresentationEvaluationContainer<VectorX<double>> container {x};
    this->evaluate<VectorX<double>>(f_unrestricted, container);

    return container.evaluation();
}


/**
 *  Calculate the matrix-vector product of (the matrix representation of) a two-electron operator with the given coefficient vector.
 *
 *  @param g                A restricted two-electron operator expressed in an orthonormal orbital basis.
 *  @param x                The coefficient vector of a linear expansion.
 *
 *  @return The coefficient vector of the linear expansion after being acted on with the given (matrix representation of) the two-electron operator.
 */
VectorX<double> SpinResolvedSelectedONVBasis::evaluateOperatorMatrixVectorProduct(const ScalarRSQTwoElectronOperator<double>& g, const VectorX<double>& x) const {

    if (g.numberOfOrbitals() != this->numberOfOrbitals()) {
        throw std::invalid_argument("SpinResolvedSelectedONVBasis::evaluateOperatorMatrixVectorProduct(const ScalarRSQTwoElectronOperator<double>&, const VectorX<double>&): The number of orbitals of this ONV basis and the operator are incompatible.");
    }

    // Use the `USQHamiltonian`'s general evaluation function, because even adding zero-valued one-electron operators won't have an impact. This would be different if we would split up the evaluation in one- and two-electron operator evaluations, which would require two times the double iterations over the whole ONV basis.
    const auto zero = ScalarUSQOneElectronOperator<double>::Zero(g.numberOfOrbitals());
    const auto g_unrestricted = ScalarUSQTwoElectronOperator<double>::FromRestricted(g);
    const USQHamiltonian<double> hamiltonian {zero, g_unrestricted};

    // Initialize a container for the matrix-vector product, and fill it with the general evaluation function.
    MatrixRepresentationEvaluationContainer<VectorX<double>> container {x};
    this->evaluate<VectorX<double>>(hamiltonian, container);

    return container.evaluation();
}


/**
 *  Calculate the matrix-vector product of (the matrix representation of) a restricted Hamiltonian with the given coefficient vector.
 *
 *  @param hamiltonian      A restricted Hamiltonian expressed in an orthonormal orbital basis.
 *  @param x                The coefficient vector of a linear expansion.
 *
 *  @return The coefficient vector of the linear expansion after being acted on with the given (matrix representation of) the Hamiltonian.
 */
VectorX<double> SpinResolvedSelectedONVBasis::evaluateOperatorMatrixVectorProduct(const RSQHamiltonian<double>& hamiltonian, const VectorX<double>& x) const {

    // By delegating the actual implementation of this method to its unrestricted counterpart, we avoid code duplication for the restricted part.
    // This does not affect performance significantly, because the bottleneck will always be the double iteration over the whole ONV basis.
    const auto h_unrestricted = ScalarUSQOneElectronOperator<double>::FromRestricted(hamiltonian.core());
    const auto g_unrestricted = ScalarUSQTwoElectronOperator<double>::FromRestricted(hamiltonian.twoElectron());
    const USQHamiltonian<double> unrestricted_hamiltonian {h_unrestricted, g_unrestricted};

    return this->evaluateOperatorMatrixVectorProduct(unrestricted_hamiltonian, x);
}


/*
 *  MARK: Sparse unrestricted operator evaluations
 */

/**
 *  Calculate the sparse matrix representation of an unrestricted Hamiltonian in this ONV basis.
 *
 *  @param hamiltonian      An unrestricted Hamiltonian expressed in an orthonormal orbital basis.
 *
 *  @return A sparse matrix represention of the Hamiltonian.
 */
Eigen::SparseMatrix<double> SpinResolvedSelectedONVBasis::evaluateOperatorSparse(const USQHamiltonian<double>& hamiltonian) const {

    if (hamiltonian.numberOfOrbitals() != this->numberOfOrbitals()) {
        throw std::invalid_argument("SpinResolvedSelectedONVBasis::evaluateOperatorSparse(const USQHamiltonian<double>&): The number of orbitals of the ONV basis and the Hamiltonian are incompatible.");
    }

    // Initialize a container for the sparse matrix representation, and reserve an appropriate amount of memory for it.
    MatrixRepresentationEvaluationContainer<Eigen::SparseMatrix<double>> container {this->dimension()};
    size_t memory = this->dimension() + this->dimension() * this->K * this->K * (this->N_alpha + this->N_beta) * (this->N_alpha + this->N_beta);
    container.reserve(memory);

    // Evaluate the Hamiltonian and add the evaluations to the sparse matrix representation.
    this->evaluate<Eigen::SparseMatrix<double>>(hamiltonian, container);

    // Finalize the creation of the sparse matrix and return the result.
    container.addToMatrix();
    return container.evaluation();
}


/*
 *  MARK: Unrestricted matrix-vector product evaluations
 */

/**
 *  Calculate the matrix-vector product of (the matrix representation of) an unrestricted Hamiltonian with the given coefficient vector.
 *
 *  @param hamiltonian      An unrestricted Hamiltonian expressed in an orthonormal orbital basis.
 *  @param x                The coefficient vector of a linear expansion.
 *
 *  @return The coefficient vector of the linear expansion after being acted on with the given (matrix representation of) the Hamiltonian.
 */
VectorX<double> SpinResolvedSelectedONVBasis::evaluateOperatorMatrixVectorProduct(const USQHamiltonian<double>& hamiltonian, const VectorX<double>& x) const {

    if (hamiltonian.numberOfOrbitals() != this->numberOfOrbitals()) {
        throw std::invalid_argument("SpinResolvedSelectedONVBasis::evaluateOperatorMatrixVectorProduct(const USQHamiltonian<double>&, const VectorX<double>& x): The number of orbitals of this ONV basis and the given Hamiltonian are incompatible.");
    }

    // Initialize a container for the matrix-vector product, and fill it with the general evaluation function.
    MatrixRepresentationEvaluationContainer<VectorX<double>> container {x};
    this->evaluate<VectorX<double>>(hamiltonian, container);

    return container.evaluation();
}


}  // namespace GQCP
