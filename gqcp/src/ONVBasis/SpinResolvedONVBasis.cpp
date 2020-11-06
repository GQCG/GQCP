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

#include "ONVBasis/SpinResolvedONVBasis.hpp"

#include <boost/math/special_functions.hpp>
#include <boost/numeric/conversion/converter.hpp>


namespace GQCP {


/*
 *  MARK: Constructors
 */

/**
 *  @param K            The number of alpha or beta spin-orbitals.
 *  @param N_alpha      The number of alpha electrons, i.e. the number of occupied alpha spin-orbitals.
 *  @param N_beta       The number of beta electrons, i.e. the number of occupied beta spin-orbitals.
 */
SpinResolvedONVBasis::SpinResolvedONVBasis(const size_t K, const size_t N_alpha, const size_t N_beta) :
    SpinResolvedBase {SpinUnresolvedONVBasis(K, N_alpha), SpinUnresolvedONVBasis(K, N_beta)} {

    // Calculate the alpha coupling elements beforehand, since this calculation is required many times in evaluating alpha-beta mixed two-electron operators.
    this->alpha_couplings = this->alpha().calculateOneElectronCouplings();
}


/*
 *  MARK: General information
 */

/**
 *  Calculate the dimension of a spin-resolved ONV basis with a given number of orbitals and electrons.
 * 
 *  @param K            The number of alpha or beta spin-orbitals.
 *  @param N_alpha      The number of alpha electrons, i.e. the number of occupied alpha spin-orbitals.
 *  @param N_beta       The number of beta electrons, i.e. the number of occupied beta spin-orbitals.
 *
 *  @return The dimension of a spin-resolved ONV basis.
 */
size_t SpinResolvedONVBasis::calculateDimension(const size_t M, const size_t N_alpha, const size_t N_beta) {

    double alpha_dim = SpinUnresolvedONVBasis::calculateDimension(M, N_alpha);
    double beta_dim = SpinUnresolvedONVBasis::calculateDimension(M, N_beta);

    try {
        return boost::numeric::converter<size_t, double>::convert(alpha_dim * beta_dim);
    } catch (boost::numeric::bad_numeric_cast& e) {
        throw std::overflow_error("SpinResolvedONVBasis::calculateDimension(size_t, size_t, size_t): " + std::string(e.what()));
    }
}


/**
 *  @return The dimension of this ONV basis.
 */
size_t SpinResolvedONVBasis::dimension() const {

    const auto K = this->alpha().numberOfOrbitals();
    const auto N_alpha = this->alpha().numberOfElectrons();
    const auto N_beta = this->beta().numberOfElectrons();

    return SpinResolvedONVBasis::calculateDimension(K, N_alpha, N_beta);
}


/*
 *  MARK: Couplings
 */

/**
 *  Calculate the one-electron operator intermediate that is required for the calculation of "theta(pq)" in Helgaker, Jørgensen, Olsen (2000). It is a partitioning of the mixed component of the unrestricted two-electron operator g(ab)_{pqrs}, resulting in a one-electron operator t(b)_{rs}.
 *
 *  @param p            The first index of the two-electron operator.
 *  @param q            The second index of the two-electron operator.
 *  @param g_ab_op      The two-electron operator.
 *
 *  @return The intermediate one-electron operator that is required for the calculation of "theta(pq)".
 */
ScalarUSQOneElectronOperatorComponent<double> SpinResolvedONVBasis::calculateOneElectronPartition(const size_t p, const size_t q, const ScalarMixedUSQTwoElectronOperatorComponent<double>& g_ab_op) const {

    // Prepare some variables.
    const auto K = g_ab_op.numberOfOrbitals();
    const auto& g_ab = g_ab_op.parameters();

    SquareMatrix<double> t = SquareMatrix<double>::Zero(K);


    // Construct the parameters of the one-electron partitioning.
    for (size_t r = 0; r < K; r++) {
        for (size_t s = 0; s < K; s++) {
            t(r, s) += g_ab(p, q, r, s);
        }
    }

    return ScalarUSQOneElectronOperatorComponent<double> {t};
}


/*
 *  MARK: Address calculations
 */

/**
 *  Calculate the compound address of an ONV represented by the two given alpha- and beta-addresses.
 * 
 *  @param I_alpha              the alpha-address
 *  @param I_beta               the beta-address
 * 
 *  @return the compound address of an ONV represented by the two given alpha- and beta-addresses.
 */
size_t SpinResolvedONVBasis::compoundAddress(const size_t I_alpha, const size_t I_beta) const {

    const auto dim_beta = this->beta().dimension();

    return I_alpha * dim_beta + I_beta;
}


/*
 *  MARK: Iterations
 */

/**
 *  Iterate over all ONVs (implicitly, by resolving in their spin components) in this ONV basis and apply the given callback function.
 * 
    *  @param callback             The function to be applied in every iteration. Its arguments are two pairs of spin-unresolved ONVs and their corresponding addresses, where the first two arguments are related to alpha-spin. The last two arguments are related to beta-spin.
 */
void SpinResolvedONVBasis::forEach(const std::function<void(const SpinUnresolvedONV&, const size_t, const SpinUnresolvedONV&, const size_t)>& callback) const {

    const auto dim_alpha = this->alpha().dimension();
    const auto dim_beta = this->beta().dimension();

    SpinUnresolvedONV onv_alpha = this->alpha().constructONVFromAddress(0);
    SpinUnresolvedONV onv_beta = this->beta().constructONVFromAddress(0);

    for (size_t Ia = 0; Ia < dim_alpha; Ia++) {  // Ia loops over addresses of alpha spin strings

        this->beta().transformONVCorrespondingToAddress(onv_beta, 0);  // reset the beta ONV to the one with the first address
        for (size_t Ib = 0; Ib < dim_beta; Ib++) {                     // Ib loops over addresses of beta spin strings

            callback(onv_alpha, Ia, onv_beta, Ib);

            if (Ib < dim_beta - 1) {  // prevent the last permutation from occurring
                this->beta().transformONVToNextPermutation(onv_beta);
            }
        }  // beta address (Ib) loop

        if (Ia < dim_alpha - 1) {  // prevent the last permutation from occurring
            this->alpha().transformONVToNextPermutation(onv_alpha);
        }
    }  // alpha address (Ia) loop
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
SquareMatrix<double> SpinResolvedONVBasis::evaluateOperatorDense(const ScalarRSQOneElectronOperator<double>& f) const {

    // Prepare some variables.
    SquareMatrix<double> F = SquareMatrix<double>::Zero(this->dimension());

    const auto dim_alpha = this->alpha().dimension();
    const auto dim_beta = this->beta().dimension();


    // The total matrix representation can be calculated from the matrix representation of the alpha- and beta-parts, but we have to place the alpha- and beta-evaluations in the correct positions in the total matrix according to the choice that alpha is 'major' and beta is 'minor'.
    const auto F_alpha = this->alpha().evaluateOperatorDense(f.alpha());
    const auto F_beta = this->beta().evaluateOperatorDense(f.beta());

    // Emplace the beta evaluations in the total matrix.
    for (size_t i = 0; i < dim_alpha; i++) {
        F.block(i * dim_beta, i * dim_beta, dim_beta, dim_beta) += F_beta;
    }

    // Emplace the alpha-evaluations in the total matrix.
    const SquareMatrix<double> ones = SquareMatrix<double>::Identity(dim_beta);
    for (size_t i = 0; i < F_alpha.cols(); i++) {
        for (size_t j = 0; j < F_alpha.cols(); j++) {
            F.block(i * dim_beta, j * dim_beta, dim_beta, dim_beta) += F_alpha(i, j) * ones;
        }
    }

    return F;
}


/**
 *  Calculate the dense matrix representation of a restricted two-electron operator in this ONV basis.
 *
 *  @param g                A restricted two-electron operator expressed in an orthonormal orbital basis.
 *
 *  @return A dense matrix represention of the two-electron operator.
 */
// SquareMatrix<double> SpinResolvedONVBasis::evaluateOperatorDense(const ScalarRSQTwoElectronOperator<double>& g) const;

/**
 *  Calculate the dense matrix representation of a restricted Hamiltonian in this ONV basis.
 *
 *  @param hamiltonian      A restricted Hamiltonian expressed in an orthonormal orbital basis.
 *
 *  @return A dense matrix represention of the Hamiltonian.
 */
// SquareMatrix<double> SpinResolvedONVBasis::evaluateOperatorDense(const RSQHamiltonian<double>& hamiltonian) const;


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
// VectorX<double> SpinResolvedONVBasis::evaluateOperatorDiagonal(const ScalarRSQOneElectronOperator<double>& f_op) const;

/**
 *  Calculate the diagonal of the matrix representation of a restricted two-electron operator in this ONV basis.
 *
 *  @param g_op             A restricted two-electron operator expressed in an orthonormal orbital basis.
 *
 *  @return The diagonal of the dense matrix represention of the two-electron operator.
 */
// VectorX<double> SpinResolvedONVBasis::evaluateOperatorDiagonal(const ScalarRSQTwoElectronOperator<double>& g_op) const;

/**
 *  Calculate the diagonal of the dense matrix representation of a restricted Hamiltonian in this ONV basis.
 *
 *  @param hamiltonian      A restricted Hamiltonian expressed in an orthonormal orbital basis.
 *
 *  @return The diagonal of the dense matrix represention of the Hamiltonian.
 */
// VectorX<double> SpinResolvedONVBasis::evaluateOperatorDiagonal(const RSQHamiltonian<double>& hamiltonian) const;


/*
 *  MARK: Restricted matrix-vector product evaluations
 */

/**
 *  Calculate the matrix-vector product of (the matrix representation of) a restricted one-electron operator with the given coefficient vector.
 *
 *  @param f                A restricted one-electron operator expressed in an orthonormal orbital basis.
 *  @param x                The coefficient vector of a linear expansion.
 *
 *  @return The coefficient vector of the linear expansion after being acted on with the given (matrix representation of) the one-electron operator.
 */
// VectorX<double> SpinResolvedONVBasis::evaluateOperatorMatrixVectorProduct(const ScalarRSQOneElectronOperator<double>& f, const VectorX<double>& x) const;

/**
 *  Calculate the matrix-vector product of (the matrix representation of) a restricted two-electron operator with the given coefficient vector.
 *
 *  @param g                A restricted two-electron operator expressed in an orthonormal orbital basis.
 *  @param x                The coefficient vector of a linear expansion.
 *
 *  @return The coefficient vector of the linear expansion after being acted on with the given (matrix representation of) the two-electron operator.
 */
// VectorX<double> SpinResolvedONVBasis::evaluateOperatorMatrixVectorProduct(const ScalarRSQTwoElectronOperator<double>& g, const VectorX<double>& x) const;

/**
 *  Calculate the matrix-vector product of (the matrix representation of) a restricted Hamiltonian with the given coefficient vector.
 *
 *  @param hamiltonian      A restricted Hamiltonian expressed in an orthonormal orbital basis.
 *  @param x                The coefficient vector of a linear expansion.
 *
 *  @return The coefficient vector of the linear expansion after being acted on with the given (matrix representation of) the Hamiltonian.
 */
// VectorX<double> SpinResolvedONVBasis::evaluateOperatorMatrixVectorProduct(const RSQHamiltonian<double>& hamiltonian, const VectorX<double>& x) const;


/*
 *  MARK: Dense unrestricted operator evaluations
 */

/**
 *  Calculate the dense matrix representation of an unrestricted Hamiltonian in this ONV basis.
 *
 *  @param hamiltonian      An unrestricted Hamiltonian expressed in an orthonormal orbital basis.
 *
 *  @return A dense matrix represention of the Hamiltonian.
 */
SquareMatrix<double> SpinResolvedONVBasis::evaluateOperatorDense(const USQHamiltonian<double>& hamiltonian) const {

    if (hamiltonian.numberOfOrbitals() != this->alpha().numberOfOrbitals()) {
        throw std::invalid_argument("SpinResolvedONVBasis::evaluateOperatorDense(const USQHamiltonian<double>&): The number of orbitals of this ONV basis and the given Hamiltonian are incompatible.");
    }


    // The total matrix representation can be calculated from the matrix representation of the alpha- and beta-parts, but we have to place the alpha- and beta-evaluations in the correct positions in the total matrix according to the choice that alpha is 'major' and beta is 'minor'.
    // In order to call the semantically correct APIs, we'll have to convert the pure alpha and pure beta part of the unrestricted Hamiltonian into a generalized representation.

    const auto& h_a = ScalarGSQOneElectronOperator<double>::FromUnrestrictedComponent(hamiltonian.core().alpha());
    const auto& g_aa = ScalarGSQTwoElectronOperator<double>::FromUnrestrictedComponent(hamiltonian.twoElectron().alphaAlpha());
    const GSQHamiltonian<double> alpha_hamiltonian {h_a, g_aa};

    const auto& h_b = ScalarGSQOneElectronOperator<double>::FromUnrestrictedComponent(hamiltonian.core().beta());
    const auto& g_bb = ScalarGSQTwoElectronOperator<double>::FromUnrestrictedComponent(hamiltonian.twoElectron().betaBeta());
    const GSQHamiltonian<double> beta_hamiltonian {h_b, g_bb};


    // Prepare some other variables.
    const auto dim_alpha = this->alpha().dimension();
    const auto dim_beta = this->beta().dimension();

    auto const& g_mixed = hamiltonian.twoElectron().alphaBeta();

    SquareMatrix<double> H = SquareMatrix<double>::Zero(this->dimension());
    const auto H_a = this->alpha().evaluateOperatorDense(alpha_hamiltonian);
    const auto H_b = this->beta().evaluateOperatorDense(beta_hamiltonian);


    // BETA separated evaluations
    for (size_t i = 0; i < dim_alpha; i++) {
        H.block(i * dim_beta, i * dim_beta, dim_beta, dim_beta) += H_b;
    }

    // ALPHA separated evaluations
    const SquareMatrix<double> ones = SquareMatrix<double>::Identity(dim_beta);
    for (int i = 0; i < H_a.cols(); i++) {
        for (int j = 0; j < H_a.cols(); j++) {
            H.block(i * dim_beta, j * dim_beta, dim_beta, dim_beta) += H_a(i, j) * ones;
        }
    }


    // MIXED evaluations
    for (size_t p = 0; p < this->alpha().numberOfOrbitals(); p++) {

        const auto& alpha_coupling = this->alphaCouplings()[p * (this->alpha().numberOfOrbitals() + this->alpha().numberOfOrbitals() + 1 - p) / 2];
        const auto& P = this->calculateOneElectronPartition(p, p, g_mixed);
        const auto& beta_two_electron_intermediate = this->beta().evaluateOperatorDense(P);

        for (int i = 0; i < alpha_coupling.outerSize(); ++i) {
            for (Eigen::SparseMatrix<double>::InnerIterator it {alpha_coupling, i}; it; ++it) {
                // it.value sigma(pp) element multiplied with the sparse matrix theta(pp) : beta_two_electron_intermediate
                H.block(it.row() * dim_beta, it.col() * dim_beta, dim_beta, dim_beta) += it.value() * beta_two_electron_intermediate;
            }
        }

        for (size_t q = p + 1; q < this->alpha().numberOfOrbitals(); q++) {

            const auto& alpha_coupling = this->alphaCouplings()[p * (this->alpha().numberOfOrbitals() + this->alpha().numberOfOrbitals() + 1 - p) / 2 + q - p];
            const auto& P = calculateOneElectronPartition(p, q, g_mixed);
            const auto& beta_two_electron_intermediate = this->beta().evaluateOperatorDense(P);

            for (int i = 0; i < alpha_coupling.outerSize(); ++i) {
                for (Eigen::SparseMatrix<double>::InnerIterator it {alpha_coupling, i}; it; ++it) {
                    // it.value (sigma(pq) + sigma(qp)) element multiplied with the sparse matrix theta(pq) : beta_two_electron_intermediate
                    H.block(it.row() * dim_beta, it.col() * dim_beta, dim_beta, dim_beta) += it.value() * beta_two_electron_intermediate;
                }
            }
        }
    }

    return H;
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
// VectorX<double> SpinResolvedONVBasis::evaluateOperatorDiagonal(const USQHamiltonian<double>& hamiltonian) const;


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
// VectorX<double> SpinResolvedONVBasis::evaluateOperatorMatrixVectorProduct(const USQHamiltonian<double>& usq_hamiltonian, const VectorX<double>& x) const;


/**
 *  Calculate the dense matrix representation of a restricted two-electron operator in this ONV basis.
 *
 *  @param g                A restricted two-electron operator expressed in an orthonormal orbital basis.
 *
 *  @return A dense matrix represention of the two-electron operator.
 */
// SquareMatrix<double> SpinResolvedONVBasis::evaluateOperatorDense(const ScalarRSQTwoElectronOperator<double>& g) const {

//     // Prepare some variables.
//     SquareMatrix<double> G = SquareMatrix<double>::Zero(this->dimension());

//     const auto dim_alpha = this->alpha().dimension();
//     const auto dim_beta = this->beta().dimension();


//     // For the two-electron evaluations, we have 'pure' combinations (i.e. alpha-alpha or beta-beta), and 'mixed' combinations (i.e. alpha-beta or beta-alpha).
//     const auto G_a = this->alpha().evaluateOperatorDense(two_op.alphaAlpha());
//     const auto G_b = this->beta().evaluateOperatorDense(two_op.betaBeta());


//     // BETA separated evaluations
//     for (size_t i = 0; i < dim_alpha; i++) {
//         G.block(i * dim_beta, i * dim_beta, dim_beta, dim_beta) += G_b;
//     }

//     // ALPHA separated evaluations
//     const SquareMatrix<double> ones = SquareMatrix<double>::Identity(dim_beta);
//     for (int i = 0; i < alpha_evaluation.cols(); i++) {
//         for (int j = 0; j < alpha_evaluation.cols(); j++) {
//             G.block(i * dim_beta, j * dim_beta, dim_beta, dim_beta) += G_a(i, j) * ones;
//         }
//     }

//     // MIXED evaluations
//     for (size_t p = 0; p < this->alpha().numberOfOrbitals(); p++) {

//         const auto& alpha_coupling = this->alphaCouplings()[p * (this->alpha().numberOfOrbitals() + this->alpha().numberOfOrbitals() + 1 - p) / 2];
//         const auto& P = this->calculateOneElectronPartition(p, p, two_op);
//         const auto& beta_two_electron_intermediate = this->beta().evaluateOperatorDense(P);

//         for (int i = 0; i < alpha_coupling.outerSize(); ++i) {
//             for (Eigen::SparseMatrix<double>::InnerIterator it {alpha_coupling, i}; it; ++it) {
//                 // it.value sigma(pp) element multiplied with the sparse matrix theta(pp) : beta_two_electron_intermediate
//                 G.block(it.row() * dim_beta, it.col() * dim_beta, dim_beta, dim_beta) += it.value() * beta_two_electron_intermediate;
//             }
//         }

//         for (size_t q = p + 1; q < this->alpha().numberOfOrbitals(); q++) {

//             const auto& alpha_coupling = this->alphaCouplings()[p * (this->alpha().numberOfOrbitals() + this->alpha().numberOfOrbitals() + 1 - p) / 2 + q - p];
//             const auto& P = this->calculateOneElectronPartition(p, q, two_op);
//             const auto& beta_two_electron_intermediate = this->beta().evaluateOperatorDense(P);

//             for (int i = 0; i < alpha_coupling.outerSize(); ++i) {
//                 for (Eigen::SparseMatrix<double>::InnerIterator it {alpha_coupling, i}; it; ++it) {
//                     // it.value (sigma(pq) + sigma(qp)) element multiplied with the sparse matrix theta(pq) : beta_two_electron_intermediate
//                     G.block(it.row() * dim_beta, it.col() * dim_beta, dim_beta, dim_beta) += it.value() * beta_two_electron_intermediate;
//                 }
//             }
//         }
//     }

//     return G;
// }


/**
 *  Evaluate the Hamiltonian in a dense matrix
 *
 *  @param sq_hamiltonian               the Hamiltonian expressed in an orthonormal basis
 *  @param diagonal_values              bool to indicate if diagonal values will be calculated
 *
 *  @return the Hamiltonian's evaluation in a dense matrix with the dimensions of the spin-resolved ONV basis
 */
// SquareMatrix<double> SpinResolvedONVBasis::evaluateOperatorDense(const RSQHamiltonian<double>& sq_hamiltonian, const bool diagonal_values) const {

//     SquareMatrix<double> total_evaluation = SquareMatrix<double>::Zero(this->dimension());

//     auto dim_alpha = this->alpha().dimension();
//     auto dim_beta = this->beta().dimension();

//     auto beta_evaluation = this->beta().evaluateOperatorDense(sq_hamiltonian, diagonal_values);
//     auto alpha_evaluation = this->alpha().evaluateOperatorDense(sq_hamiltonian, diagonal_values);

//     // BETA separated evaluations
//     for (size_t i = 0; i < dim_alpha; i++) {
//         total_evaluation.block(i * dim_beta, i * dim_beta, dim_beta, dim_beta) += beta_evaluation;
//     }

//     // ALPHA separated evaluations
//     const SquareMatrix<double> ones = SquareMatrix<double>::Identity(dim_beta);
//     for (int i = 0; i < alpha_evaluation.cols(); i++) {
//         for (int j = 0; j < alpha_evaluation.cols(); j++) {
//             total_evaluation.block(i * dim_beta, j * dim_beta, dim_beta, dim_beta) += alpha_evaluation(i, j) * ones;
//         }
//     }

//     // MIXED evaluations
//     for (size_t p = 0; p < this->alpha().numberOfOrbitals(); p++) {

//         const auto& alpha_coupling = this->alphaCouplings()[p * (this->alpha().numberOfOrbitals() + this->alpha().numberOfOrbitals() + 1 - p) / 2];
//         const auto& P = this->calculateOneElectronPartition(p, p, sq_hamiltonian.twoElectron());
//         const auto& beta_two_electron_intermediate = this->beta().evaluateOperatorDense(P, diagonal_values);

//         for (int i = 0; i < alpha_coupling.outerSize(); ++i) {
//             for (Eigen::SparseMatrix<double>::InnerIterator it {alpha_coupling, i}; it; ++it) {
//                 // it.value sigma(pp) element multiplied with the sparse matrix theta(pp) : beta_two_electron_intermediate
//                 total_evaluation.block(it.row() * dim_beta, it.col() * dim_beta, dim_beta, dim_beta) += it.value() * beta_two_electron_intermediate;
//             }
//         }

//         for (size_t q = p + 1; q < this->alpha().numberOfOrbitals(); q++) {

//             const auto& alpha_coupling = this->alphaCouplings()[p * (this->alpha().numberOfOrbitals() + this->alpha().numberOfOrbitals() + 1 - p) / 2 + q - p];
//             const auto& P = calculateOneElectronPartition(p, q, sq_hamiltonian.twoElectron());
//             const auto& beta_two_electron_intermediate = this->beta().evaluateOperatorDense(P, true);

//             for (int i = 0; i < alpha_coupling.outerSize(); ++i) {
//                 for (Eigen::SparseMatrix<double>::InnerIterator it {alpha_coupling, i}; it; ++it) {
//                     // it.value (sigma(pq) + sigma(qp)) element multiplied with the sparse matrix theta(pq) : beta_two_electron_intermediate
//                     total_evaluation.block(it.row() * dim_beta, it.col() * dim_beta, dim_beta, dim_beta) += it.value() * beta_two_electron_intermediate;
//                 }
//             }
//         }
//     }

//     return total_evaluation;
// }


// /**
//  *  Evaluate the diagonal of the operator in this spin-resolved ONV basis
//  *
//  *  @param one_op               the one-electron operator in an orthonormal orbital basis to be evaluated in the spin-resolved ONV basis
//  *
//  *  @return the operator's diagonal evaluation in a vector with the dimension of the spin-resolved ONV basis
//  */
// VectorX<double> SpinResolvedONVBasis::evaluateOperatorDiagonal(const ScalarRSQOneElectronOperator<double>& one_op) const {

//     const auto M = one_op.numberOfOrbitals();
//     if (M != this->alpha().numberOfOrbitals()) {
//         throw std::invalid_argument("SpinResolvedONVBasis::evaluateOperatorDiagonal(ScalarRSQOneElectronOperator<double>): Basis functions of the spin-resolved ONV basis and the operator are incompatible.");
//     }

//     const auto dim_alpha = this->alpha().dimension();
//     const auto dim_beta = this->beta().dimension();
//     const auto& one_op_par = one_op.parameters();

//     VectorX<double> diagonal = VectorX<double>::Zero(this->dim);

//     SpinUnresolvedONV onv_alpha = this->alpha().constructONVFromAddress(0);
//     SpinUnresolvedONV onv_beta = this->beta().constructONVFromAddress(0);
//     for (size_t Ia = 0; Ia < dim_alpha; Ia++) {  // Ia loops over addresses of alpha spin strings

//         this->beta().transformONVCorrespondingToAddress(onv_beta, 0);

//         for (size_t Ib = 0; Ib < dim_beta; Ib++) {  // Ib loops over addresses of beta spin strings

//             for (size_t e_a = 0; e_a < this->alpha().numberOfElectrons(); e_a++) {  // loop over alpha electrons

//                 size_t p = onv_alpha.occupationIndexOf(e_a);
//                 diagonal(Ia * dim_beta + Ib) += one_op_par(p, p);

//             }  // e_a loop

//             for (size_t e_b = 0; e_b < this->beta().numberOfElectrons(); e_b++) {  // loop over beta electrons

//                 size_t p = onv_beta.occupationIndexOf(e_b);
//                 diagonal(Ia * dim_beta + Ib) += one_op_par(p, p);
//             }

//             if (Ib < dim_beta - 1) {  // prevent the last permutation from occurring
//                 this->beta().transformONVToNextPermutation(onv_beta);
//             }
//         }  // beta address (Ib) loop

//         if (Ia < dim_alpha - 1) {  // prevent the last permutation from occurring
//             this->alpha().transformONVToNextPermutation(onv_alpha);
//         }
//     }  // alpha address (Ia) loop

//     return diagonal;
// }


// /**
//  *  Evaluate the diagonal of the operator in this spin-resolved ONV basis
//  *
//  *  @param two_op               the two-electron operator in an orthonormal orbital basis to be evaluated in the spin-resolved ONV basis
//  *
//  *  @return the operator's diagonal evaluation in a vector with the dimension of the spin-resolved ONV basis
//  */
// VectorX<double> SpinResolvedONVBasis::evaluateOperatorDiagonal(const ScalarRSQTwoElectronOperator<double>& two_op) const {

//     const auto M = two_op.numberOfOrbitals();
//     if (M != this->alpha().numberOfOrbitals()) {
//         throw std::invalid_argument("SpinResolvedONVBasis::evaluateOperatorDiagonal(ScalarRSQTwoElectronOperator<double>): Basis functions of the SpinUnresolvedONV basis and the operator are incompatible.");
//     }

//     const auto dim_alpha = this->alpha().dimension();
//     const auto dim_beta = this->beta().dimension();
//     const auto& two_op_par = two_op.parameters();
//     const auto k = two_op.effectiveOneElectronPartition().parameters();

//     // Diagonal contributions
//     VectorX<double> diagonal = VectorX<double>::Zero(this->dim);

//     SpinUnresolvedONV onv_alpha = this->alpha().constructONVFromAddress(0);
//     SpinUnresolvedONV onv_beta = this->beta().constructONVFromAddress(0);
//     for (size_t Ia = 0; Ia < dim_alpha; Ia++) {  // Ia loops over addresses of alpha spin strings

//         this->beta().transformONVCorrespondingToAddress(onv_beta, 0);

//         for (size_t Ib = 0; Ib < dim_beta; Ib++) {  // Ib loops over addresses of beta spin strings

//             for (size_t e_a = 0; e_a < this->alpha().numberOfElectrons(); e_a++) {  // loop over alpha electrons

//                 size_t p = onv_alpha.occupationIndexOf(e_a);
//                 diagonal(Ia * dim_beta + Ib) += k(p, p);

//                 for (size_t q = 0; q < this->alpha().numberOfOrbitals(); q++) {  // q loops over SOs
//                     if (onv_alpha.isOccupied(q)) {      // q is in Ia
//                         diagonal(Ia * dim_beta + Ib) += 0.5 * two_op_par(p, p, q, q);
//                     } else {  // q is not in I_alpha
//                         diagonal(Ia * dim_beta + Ib) += 0.5 * two_op_par(p, q, q, p);
//                     }

//                     if (onv_beta.isOccupied(q)) {  // q is in Ib
//                         diagonal(Ia * dim_beta + Ib) += two_op_par(p, p, q, q);
//                     }
//                 }  // q loop
//             }      // e_a loop

//             for (size_t e_b = 0; e_b < this->beta().numberOfElectrons(); e_b++) {  // loop over beta electrons

//                 size_t p = onv_beta.occupationIndexOf(e_b);
//                 diagonal(Ia * dim_beta + Ib) += k(p, p);

//                 for (size_t q = 0; q < this->alpha().numberOfOrbitals(); q++) {  // q loops over SOs
//                     if (onv_beta.isOccupied(q)) {       // q is in Ib
//                         diagonal(Ia * dim_beta + Ib) += 0.5 * two_op_par(p, p, q, q);

//                     } else {  // q is not in I_beta
//                         diagonal(Ia * dim_beta + Ib) += 0.5 * two_op_par(p, q, q, p);
//                     }
//                 }  // q loop
//             }      // e_b loop

//             if (Ib < dim_beta - 1) {  // prevent the last permutation from occurring
//                 this->beta().transformONVToNextPermutation(onv_beta);
//             }
//         }  // beta address (Ib) loop

//         if (Ia < dim_alpha - 1) {  // prevent the last permutation from occurring
//             this->alpha().transformONVToNextPermutation(onv_alpha);
//         }
//     }  // alpha address (Ia) loop

//     return diagonal;
// }


// /**
//  *  Evaluate the diagonal of the Hamiltonian in this spin-resolved ONV basis
//  *
//  *  @param sq_hamiltonian           the Hamiltonian expressed in an orthonormal basis
//  *
//  *  @return the Hamiltonian's diagonal evaluation in a vector with the dimension of the spin-resolved ONV basis
//  */
// VectorX<double> SpinResolvedONVBasis::evaluateOperatorDiagonal(const RSQHamiltonian<double>& sq_hamiltonian) const {
//     return this->evaluateOperatorDiagonal(sq_hamiltonian.core()) + this->evaluateOperatorDiagonal(sq_hamiltonian.twoElectron());
// }


// /**
//  *  Evaluate the diagonal of the Hamiltonian
//  *
//  *  @param usq_hamiltonian                the Hamiltonian expressed in an unrestricted orthonormal basis
//  *
//  *  @return the Hamiltonian's diagonal evaluation in a vector with the dimension of the spin-resolved ONV basis
//  */
// VectorX<double> SpinResolvedONVBasis::evaluateOperatorDiagonal(const USQHamiltonian<double>& usq_hamiltonian) const {

//     const auto M = usq_hamiltonian.numberOfOrbitals() / 2;

//     if (!usq_hamiltonian.areSpinHamiltoniansOfSameDimension()) {
//         throw std::invalid_argument("SpinResolvedONVBasis::evaluateOperatorDiagonal(USQHamiltonian<double>): Underlying spin Hamiltonians are not of the same dimension, and this is currently required for this method");
//     }

//     if (M != this->alpha().numberOfOrbitals()) {
//         throw std::invalid_argument("SpinResolvedONVBasis::evaluateOperatorDiagonal(USQHamiltonian<double>): Basis functions of the spin-resolved ONV basis and the operator are incompatible.");
//     }

//     // Evaluation environment
//     auto const& alpha_hamiltonian = usq_hamiltonian.spinHamiltonian(Spin::alpha);
//     auto const& beta_hamiltonian = usq_hamiltonian.spinHamiltonian(Spin::beta);
//     auto const& g_mixed = usq_hamiltonian.twoElectronMixed();

//     const auto dim_alpha = this->alpha().dimension();
//     const auto dim_beta = this->beta().dimension();
//     auto k_alpha = alpha_hamiltonian.core().parameters();
//     auto k_beta = beta_hamiltonian.core().parameters();
//     const auto& two_op_par_alpha = alpha_hamiltonian.twoElectron().parameters();
//     const auto& two_op_par_beta = beta_hamiltonian.twoElectron().parameters();

//     k_alpha = k_alpha + alpha_hamiltonian.twoElectron().effectiveOneElectronPartition().parameters();
//     k_beta = k_beta + beta_hamiltonian.twoElectron().effectiveOneElectronPartition().parameters();

//     // The two_op_par_mixed variable stored as g_aabb, for integrals derived from g_bbaa we reverse the indices as follows : g_aabb(pqrs) = g_bbaa(rspq)
//     const auto& two_op_par_mixed = g_mixed.parameters();

//     VectorX<double> diagonal = VectorX<double>::Zero(this->dim);

//     SpinUnresolvedONV onv_alpha = this->alpha().constructONVFromAddress(0);
//     SpinUnresolvedONV onv_beta = this->beta().constructONVFromAddress(0);
//     for (size_t Ia = 0; Ia < dim_alpha; Ia++) {  // Ia loops over addresses of alpha spin strings

//         this->beta().transformONVCorrespondingToAddress(onv_beta, 0);

//         for (size_t Ib = 0; Ib < dim_beta; Ib++) {  // Ib loops over addresses of beta spin strings

//             for (size_t e_a = 0; e_a < this->alpha().numberOfElectrons(); e_a++) {  // loop over alpha electrons

//                 size_t p = onv_alpha.occupationIndexOf(e_a);
//                 diagonal(Ia * dim_beta + Ib) += k_alpha(p, p);

//                 for (size_t q = 0; q < this->alpha().numberOfOrbitals(); q++) {  // q loops over SOs
//                     if (onv_alpha.isOccupied(q)) {      // q is in Ia
//                         diagonal(Ia * dim_beta + Ib) += 0.5 * two_op_par_alpha(p, p, q, q);
//                     } else {  // q is not in I_alpha
//                         diagonal(Ia * dim_beta + Ib) += 0.5 * two_op_par_alpha(p, q, q, p);
//                     }

//                     if (onv_beta.isOccupied(q)) {  // q is in Ib
//                         diagonal(Ia * dim_beta + Ib) += two_op_par_mixed(p, p, q, q);
//                     }
//                 }  // q loop
//             }      // e_a loop

//             for (size_t e_b = 0; e_b < this->beta().numberOfElectrons(); e_b++) {  // loop over beta electrons

//                 size_t p = onv_beta.occupationIndexOf(e_b);
//                 diagonal(Ia * dim_beta + Ib) += k_beta(p, p);

//                 for (size_t q = 0; q < this->alpha().numberOfOrbitals(); q++) {  // q loops over SOs
//                     if (onv_beta.isOccupied(q)) {       // q is in Ib
//                         diagonal(Ia * dim_beta + Ib) += 0.5 * two_op_par_beta(p, p, q, q);

//                     } else {  // q is not in I_beta
//                         diagonal(Ia * dim_beta + Ib) += 0.5 * two_op_par_beta(p, q, q, p);
//                     }
//                 }  // q loop
//             }      // e_b loop

//             if (Ib < dim_beta - 1) {  // prevent the last permutation from occurring
//                 this->beta().transformONVToNextPermutation(onv_beta);
//             }
//         }  // beta address (Ib) loop

//         if (Ia < dim_alpha - 1) {  // prevent the last permutation from occurring
//             this->alpha().transformONVToNextPermutation(onv_alpha);
//         }
//     }  // alpha address (Ia) loop

//     return diagonal;
// }


// /**
//  *  Evaluate a one electron operator in a matrix vector product
//  *
//  *  @param one_op                       the one electron operator expressed in an orthonormal basis
//  *  @param x                            the vector upon which the evaluation acts
//  *  @param diagonal                     the diagonal evaluated in the spin-resolved ONV basis
//  *
//  *  @return the one electron operator's matrix vector product in a vector with the dimensions of the spin-resolved ONV basis
//  */
// VectorX<double> SpinResolvedONVBasis::evaluateOperatorMatrixVectorProduct(const ScalarRSQOneElectronOperator<double>& one_op, const VectorX<double>& x, const VectorX<double>& diagonal) const {

//     auto M = one_op.numberOfOrbitals();
//     if (M != this->alpha().numberOfOrbitals()) {
//         throw std::invalid_argument("SpinResolvedONVBasis::evaluateOperatorMatrixVectorProduct(ScalarRSQOneElectronOperator<double>, VectorX<double>, VectorX<double>): Basis functions of the spin-resolved ONV basis and the operator are incompatible.");
//     }

//     // Environment for evaluations
//     SpinUnresolvedONVBasis this->alpha() = this->alpha();
//     SpinUnresolvedONVBasis this->beta() = this->beta

//     const auto& alpha_couplings = this->alphaCouplings();

//     auto dim_alpha = this->alpha().dimension();
//     auto dim_beta = this->beta().dimension();

//     VectorX<double> matvec = diagonal.cwiseProduct(x);

//     // Map vector to matrix for vectorized multiplications
//     Eigen::Map<Eigen::MatrixXd> matvecmap(matvec.data(), dim_beta, dim_alpha);
//     Eigen::Map<const Eigen::MatrixXd> xmap(x.data(), dim_beta, dim_alpha);

//     // Spin-resolved evaluation
//     auto beta_evaluation = this->beta().evaluateOperatorSparse(one_op, false);
//     auto alpha_evaluation = this->alpha().evaluateOperatorSparse(one_op, false);

//     // Perform the "matvec"
//     matvecmap += xmap * alpha_evaluation + beta_evaluation * xmap;

//     return matvec;
// }


// /**
//  *  Evaluate a two electron operator in a matrix vector product
//  *
//  *  @param two_op                       the two electron operator expressed in an orthonormal basis
//  *  @param x                            the vector upon which the evaluation acts
//  *  @param diagonal                     the diagonal evaluated in the spin-resolved ONV basis
//  *
//  *  @return the two electron operator's matrix vector product in a vector with the dimensions of the spin-resolved ONV basis
//  */
// VectorX<double> SpinResolvedONVBasis::evaluateOperatorMatrixVectorProduct(const ScalarRSQTwoElectronOperator<double>& two_op, const VectorX<double>& x, const VectorX<double>& diagonal) const {

//     auto M = two_op.numberOfOrbitals();
//     if (M != this->alpha().numberOfOrbitals()) {
//         throw std::invalid_argument("SpinResolvedONVBasis::evaluateOperatorMatrixVectorProduct(ScalarRSQTwoElectronOperator<double>, VectorX<double>, VectorX<double>): Basis functions of the spin-resolved ONV basis and the operator are incompatible.");
//     }

//     // Environment for evaluations
//     SpinUnresolvedONVBasis this->alpha() = this->alpha();
//     SpinUnresolvedONVBasis this->beta() = this->beta

//     const auto& alpha_couplings = this->alphaCouplings();

//     auto dim_alpha = this->alpha().dimension();
//     auto dim_beta = this->beta().dimension();

//     VectorX<double> matvec = diagonal.cwiseProduct(x);

//     Eigen::Map<Eigen::MatrixXd> matvecmap(matvec.data(), dim_beta, dim_alpha);
//     Eigen::Map<const Eigen::MatrixXd> xmap(x.data(), dim_beta, dim_alpha);

//     // Mixed-spin evaluation
//     for (size_t p = 0; p < this->alpha().numberOfOrbitals(); p++) {

//         const auto& P = this->calculateOneElectronPartition(p, p, two_op);
//         const auto& beta_two_electron_intermediate = this->beta().evaluateOperatorSparse(P, false);

//         // matvec : sigma(pp) * X * theta(pp)
//         matvecmap += beta_two_electron_intermediate * (xmap * alpha_couplings[p * (this->alpha().numberOfOrbitals() + this->alpha().numberOfOrbitals() + 1 - p) / 2]);
//         for (size_t q = p + 1; q < this->alpha().numberOfOrbitals(); q++) {

//             const auto& P = this->calculateOneElectronPartition(p, q, two_op);
//             const auto& beta_two_electron_intermediate = this->beta().evaluateOperatorSparse(P, true);

//             // matvec : (sigma(pq) + sigma(qp)) * X * theta(pq)
//             matvecmap += beta_two_electron_intermediate * (xmap * alpha_couplings[p * (this->alpha().numberOfOrbitals() + this->alpha().numberOfOrbitals() + 1 - p) / 2 + q - p]);
//         }
//     }

//     // Spin-resolved evaluation
//     auto beta_evaluation = this->beta().evaluateOperatorSparse(two_op, false);
//     auto alpha_evaluation = this->alpha().evaluateOperatorSparse(two_op, false);

//     matvecmap += beta_evaluation * xmap + xmap * alpha_evaluation;

//     return matvec;
// }


// /**
//  *  Evaluate the Hamiltonian in a matrix vector product
//  *
//  *  @param sq_hamiltonian               the Hamiltonian expressed in an orthonormal basis
//  *  @param x                            the vector upon which the evaluation acts
//  *  @param diagonal                     the diagonal evaluated in the spin-resolved ONV basis
//  *
//  *  @return the Hamiltonian's matrix vector product in a vector with the dimensions of the spin-resolved ONV basis
//  */
// VectorX<double> SpinResolvedONVBasis::evaluateOperatorMatrixVectorProduct(const RSQHamiltonian<double>& sq_hamiltonian, const VectorX<double>& x, const VectorX<double>& diagonal) const {

//     auto M = sq_hamiltonian.numberOfOrbitals();
//     if (M != this->alpha().numberOfOrbitals()) {
//         throw std::invalid_argument("SpinResolvedONVBasis::evaluateOperatorMatrixVectorProduct(RSQHamiltonian<double>, VectorX<double>, VectorX<double>): Basis functions of the spin-resolved ONV basis and the operator are incompatible.");
//     }

//     // Environment for evaluations
//     const SpinUnresolvedONVBasis& this->alpha() = this->alpha();
//     const SpinUnresolvedONVBasis& this->beta() = this->beta

//     const auto& alpha_couplings = this->alphaCouplings();

//     auto dim_alpha = this->alpha().dimension();
//     auto dim_beta = this->beta().dimension();

//     VectorX<double> matvec = diagonal.cwiseProduct(x);

//     Eigen::Map<Eigen::MatrixXd> matvecmap(matvec.data(), dim_beta, dim_alpha);
//     Eigen::Map<const Eigen::MatrixXd> xmap(x.data(), dim_beta, dim_alpha);

//     // Mixed-spin evaluation
//     for (size_t p = 0; p < this->alpha().numberOfOrbitals(); p++) {

//         const auto& P = this->calculateOneElectronPartition(p, p, sq_hamiltonian.twoElectron());
//         const auto& beta_two_electron_intermediate = this->beta().evaluateOperatorSparse(P, false);

//         // matvec : sigma(pp) * X * theta(pp)
//         matvecmap += beta_two_electron_intermediate * (xmap * alpha_couplings[p * (this->alpha().numberOfOrbitals() + this->alpha().numberOfOrbitals() + 1 - p) / 2]);
//         for (size_t q = p + 1; q < this->alpha().numberOfOrbitals(); q++) {

//             const auto& P = this->calculateOneElectronPartition(p, q, sq_hamiltonian.twoElectron());
//             const auto& beta_two_electron_intermediate = this->beta().evaluateOperatorSparse(P, true);

//             // matvec : (sigma(pq) + sigma(qp)) * X * theta(pq)
//             matvecmap += beta_two_electron_intermediate * (xmap * alpha_couplings[p * (this->alpha().numberOfOrbitals() + this->alpha().numberOfOrbitals() + 1 - p) / 2 + q - p]);
//         }
//     }

//     // Spin-resolved evaluation
//     auto beta_hamiltonian = this->beta().evaluateOperatorSparse(sq_hamiltonian, false);
//     auto alpha_hamiltonian = this->alpha().evaluateOperatorSparse(sq_hamiltonian, false);

//     matvecmap += beta_hamiltonian * xmap + xmap * alpha_hamiltonian;

//     return matvec;
// }


// /**
//  *  Evaluate the unrestricted Hamiltonian in a matrix vector product
//  *
//  *  @param usq_hamiltonian                the Hamiltonian expressed in an unrestricted orthonormal basis
//  *  @param x                              the vector upon which the evaluation acts
//  *  @param diagonal                       the diagonal evaluated in the spin-resolved ONV basis
//  *
//  *  @return the Hamiltonian's evaluation in a dense matrix with the dimensions of the spin-resolved ONV basis
//  */
// VectorX<double> SpinResolvedONVBasis::evaluateOperatorMatrixVectorProduct(const USQHamiltonian<double>& usq_hamiltonian, const VectorX<double>& x, const VectorX<double>& diagonal) const {

//     auto M = usq_hamiltonian.numberOfOrbitals() / 2;

//     if (!usq_hamiltonian.areSpinHamiltoniansOfSameDimension()) {
//         throw std::invalid_argument("SpinResolvedONVBasis::evaluateOperatorMatrixVectorProduct(USQHamiltonian<double>, VectorX<double> , VectorX<double>): Underlying spin Hamiltonians are not of the same dimension, and this is currently required for this method");
//     }

//     if (M != this->numberOfOrbitals()) {
//         throw std::invalid_argument("SpinResolvedONVBasis::evaluateOperatorMatrixVectorProduct(USQHamiltonian<double>, VectorX<double>, VectorX<double>): Basis functions of the spin-resolved ONV basis and usq_hamiltonian are incompatible.");
//     }

//     // Environment for evaluations
//     const SpinUnresolvedONVBasis& this->alpha() = this->alpha();
//     const SpinUnresolvedONVBasis& this->beta() = this->beta

//     const auto& alpha_couplings = this->alphaCouplings();

//     auto dim_alpha = this->alpha().dimension();
//     auto dim_beta = this->beta().dimension();

//     VectorX<double> matvec = diagonal.cwiseProduct(x);

//     Eigen::Map<Eigen::MatrixXd> matvecmap(matvec.data(), dim_beta, dim_alpha);
//     Eigen::Map<const Eigen::MatrixXd> xmap(x.data(), dim_beta, dim_alpha);

//     for (size_t p = 0; p < this->alpha().numberOfOrbitals(); p++) {

//         const auto& P = this->calculateOneElectronPartition(p, p, usq_hamiltonian.twoElectronMixed());
//         const auto& beta_two_electron_intermediate = this->beta().evaluateOperatorDense(P, false);

//         // sigma(pp) * X * theta(pp)
//         matvecmap += beta_two_electron_intermediate * (xmap * alpha_couplings[p * (this->alpha().numberOfOrbitals() + this->alpha().numberOfOrbitals() + 1 - p) / 2]);
//         for (size_t q = p + 1; q < this->alpha().numberOfOrbitals(); q++) {

//             const auto& P = this->calculateOneElectronPartition(p, q, usq_hamiltonian.twoElectronMixed());
//             const auto& beta_two_electron_intermediate = this->beta().evaluateOperatorDense(P, true);

//             // (sigma(pq) + sigma(qp)) * X * theta(pq)
//             matvecmap += beta_two_electron_intermediate * (xmap * alpha_couplings[p * (this->alpha().numberOfOrbitals() + this->alpha().numberOfOrbitals() + 1 - p) / 2 + q - p]);
//         }
//     }

//     auto beta_hamiltonian = this->beta().evaluateOperatorSparse(usq_hamiltonian.spinHamiltonian(Spin::beta), false);
//     auto alpha_hamiltonian = this->alpha().evaluateOperatorSparse(usq_hamiltonian.spinHamiltonian(Spin::alpha), false);

//     matvecmap += beta_hamiltonian * xmap + xmap * alpha_hamiltonian;

//     return matvec;
// }


}  // namespace GQCP
