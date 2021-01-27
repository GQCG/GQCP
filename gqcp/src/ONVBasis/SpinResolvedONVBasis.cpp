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

    for (size_t Ia = 0; Ia < dim_alpha; Ia++) {  // Ia loops over addresses of alpha ONVs

        this->beta().transformONVCorrespondingToAddress(onv_beta, 0);  // reset the beta ONV to the one with the first address
        for (size_t Ib = 0; Ib < dim_beta; Ib++) {                     // Ib loops over addresses of beta ONVs

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
SquareMatrix<double> SpinResolvedONVBasis::evaluateOperatorDense(const ScalarRSQTwoElectronOperator<double>& g) const {

    // We choose to avoid code duplication by evaluating an equivalent restricted Hamiltonian with zero-valued core contributions.
    const auto zero = ScalarRSQOneElectronOperator<double>::Zero(g.numberOfOrbitals());
    const RSQHamiltonian<double> hamiltonian {zero, g};

    return this->evaluateOperatorDense(hamiltonian);
}


/**
 *  Calculate the dense matrix representation of a restricted Hamiltonian in this ONV basis.
 *
 *  @param hamiltonian      A restricted Hamiltonian expressed in an orthonormal orbital basis.
 *
 *  @return A dense matrix represention of the Hamiltonian.
 */
SquareMatrix<double> SpinResolvedONVBasis::evaluateOperatorDense(const RSQHamiltonian<double>& hamiltonian) const {

    // We can avoid code duplication by delegating this method to the evaluation of an unrestricted Hamiltonian. This entails a small speed decrease, but doesn't change the order of the scaling.
    const auto h_unrestricted = ScalarUSQOneElectronOperator<double>::FromRestricted(hamiltonian.core());
    const auto g_unrestricted = ScalarUSQTwoElectronOperator<double>::FromRestricted(hamiltonian.twoElectron());
    const USQHamiltonian<double> unrestricted_hamiltonian {h_unrestricted, g_unrestricted};

    return this->evaluateOperatorDense(unrestricted_hamiltonian);
}


/**
 *  Calculate the dense matrix representation of a Hubbard Hamiltonian in this ONV basis.
 *
 *  @param hamiltonian      A Hubbard Hamiltonian expressed in an orthonormal orbital basis.
 *
 *  @return A dense matrix represention of the Hamiltonian.
 */
SquareMatrix<double> SpinResolvedONVBasis::evaluateOperatorDense(const HubbardHamiltonian<double>& hamiltonian) const {

    if (hamiltonian.numberOfLatticeSites() != this->alpha().numberOfOrbitals()) {
        throw std::invalid_argument("SpinResolvedONVBasis::evaluateOperatorDense(const HubbardHamiltonian<double>&): The number of spatial orbitals of this ONV basis and the number of lattice sites for the Hubbard Hamiltonian are incompatible.");
    }


    // Split up the evaluation of the Hamiltonian in one- and two-electron contributions.
    // The hopping terms lead to one-electron contributions.
    // The two-electron on-site repulsion terms only lead to diagonal contributions.
    const auto H1 = this->evaluateOperatorDense(hamiltonian.core());
    const SquareMatrix<double> H2 {this->evaluateOperatorDiagonal(hamiltonian).asDiagonal()};

    return H1 + H2;
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
VectorX<double> SpinResolvedONVBasis::evaluateOperatorDiagonal(const ScalarRSQOneElectronOperator<double>& f_op) const {

    if (f_op.numberOfOrbitals() != this->alpha().numberOfOrbitals()) {
        throw std::invalid_argument("SpinResolvedONVBasis::evaluateOperatorDiagonal(ScalarRSQOneElectronOperator<double>): The number of orbitals of the spin-resolved ONV basis and the operator are incompatible.");
    }

    // Prepare some variables.
    const auto dim_alpha = this->alpha().dimension();
    const auto dim_beta = this->beta().dimension();
    const auto& f = f_op.parameters();

    VectorX<double> diagonal = VectorX<double>::Zero(this->dimension());

    SpinUnresolvedONV onv_alpha = this->alpha().constructONVFromAddress(0);
    SpinUnresolvedONV onv_beta = this->beta().constructONVFromAddress(0);
    for (size_t Ia = 0; Ia < dim_alpha; Ia++) {  // Ia loops over addresses of alpha ONVs

        this->beta().transformONVCorrespondingToAddress(onv_beta, 0);
        for (size_t Ib = 0; Ib < dim_beta; Ib++) {  // Ib loops over addresses of beta ONVs

            for (size_t e_a = 0; e_a < this->alpha().numberOfElectrons(); e_a++) {  // loop over alpha electrons
                size_t p = onv_alpha.occupationIndexOf(e_a);
                diagonal(Ia * dim_beta + Ib) += f(p, p);
            }  // e_a loop

            for (size_t e_b = 0; e_b < this->beta().numberOfElectrons(); e_b++) {  // loop over beta electrons
                size_t p = onv_beta.occupationIndexOf(e_b);
                diagonal(Ia * dim_beta + Ib) += f(p, p);
            }

            if (Ib < dim_beta - 1) {  // prevent the last permutation from occurring
                this->beta().transformONVToNextPermutation(onv_beta);
            }
        }  // beta address (Ib) loop

        if (Ia < dim_alpha - 1) {  // prevent the last permutation from occurring
            this->alpha().transformONVToNextPermutation(onv_alpha);
        }
    }  // alpha address (Ia) loop

    return diagonal;
}


/**
 *  Calculate the diagonal of the matrix representation of a restricted two-electron operator in this ONV basis.
 *
 *  @param g                A restricted two-electron operator expressed in an orthonormal orbital basis.
 *
 *  @return The diagonal of the dense matrix represention of the two-electron operator.
 */
VectorX<double> SpinResolvedONVBasis::evaluateOperatorDiagonal(const ScalarRSQTwoElectronOperator<double>& g) const {

    // We choose to avoid code duplication by evaluating an equivalent restricted Hamiltonian with zero-valued core contributions.
    const auto zero = ScalarRSQOneElectronOperator<double>::Zero(g.numberOfOrbitals());
    const RSQHamiltonian<double> hamiltonian {zero, g};

    return this->evaluateOperatorDiagonal(hamiltonian);
}


/**
 *  Calculate the diagonal of the dense matrix representation of a restricted Hamiltonian in this ONV basis.
 *
 *  @param hamiltonian      A restricted Hamiltonian expressed in an orthonormal orbital basis.
 *
 *  @return The diagonal of the dense matrix represention of the Hamiltonian.
 */
VectorX<double> SpinResolvedONVBasis::evaluateOperatorDiagonal(const RSQHamiltonian<double>& hamiltonian) const {

    // We can avoid code duplication by delegating this method to the evaluation of an unrestricted Hamiltonian. This entails a small speed decrease, but doesn't change the order of the scaling.
    const auto h_unrestricted = ScalarUSQOneElectronOperator<double>::FromRestricted(hamiltonian.core());
    const auto g_unrestricted = ScalarUSQTwoElectronOperator<double>::FromRestricted(hamiltonian.twoElectron());
    const USQHamiltonian<double> unrestricted_hamiltonian {h_unrestricted, g_unrestricted};

    return this->evaluateOperatorDiagonal(unrestricted_hamiltonian);
}


/**
 *  Calculate the diagonal of the dense matrix representation of a Hubbard Hamiltonian in this ONV basis.
 *
 *  @param hamiltonian      A Hubbard Hamiltonian expressed in an orthonormal orbital basis.
 *
 *  @return The diagonal of the dense matrix represention of the Hamiltonian.
 */
VectorX<double> SpinResolvedONVBasis::evaluateOperatorDiagonal(const HubbardHamiltonian<double>& hamiltonian) const {

    if (hamiltonian.numberOfLatticeSites() != this->alpha().numberOfOrbitals()) {
        throw std::invalid_argument("SpinResolvedONVBasis::evaluateOperatorDiagonal(const HubbardHamiltonian<double>&): The number of spatial orbitals of this ONV basis and the number of lattice sites for the Hubbard Hamiltonian are incompatible.");
    }


    // Prepare some variables.
    const auto dim_alpha = this->alpha().dimension();
    const auto dim_beta = this->beta().dimension();
    const auto& H = hamiltonian.hoppingMatrix().matrix();


    // Calculate the diagonal contributions resulting from the two-electron on-site interactions by iterating over all ONVs.
    VectorX<double> diagonal = VectorX<double>::Zero(this->dimension());

    SpinUnresolvedONV onv_alpha = this->alpha().constructONVFromAddress(0);
    SpinUnresolvedONV onv_beta = this->beta().constructONVFromAddress(0);
    for (size_t Ia = 0; Ia < dim_alpha; Ia++) {  // Ia loops over the addresses of alpha ONVs.
        this->beta().transformONVCorrespondingToAddress(onv_beta, 0);

        for (size_t Ib = 0; Ib < dim_beta; Ib++) {  // Ib loops over addresses of beta ONVs.
            const auto I = this->compoundAddress(Ia, Ib);

            // There is a contribution for all orbital indices p that are occupied both in the alpha- and beta ONV.
            std::vector<size_t> occupations = onv_alpha.findMatchingOccupations(onv_beta);
            for (const auto& p : occupations) {
                diagonal(I) += H(p, p);  // The two-electron (on-site repulsion) contributions are on the diagonal of the hopping matrix.
            }

            if (Ib < dim_beta - 1) {  // Prevent the last permutation from occurring.
                this->beta().transformONVToNextPermutation(onv_beta);
            }
        }  // Beta address (Ib) loop.

        if (Ia < dim_alpha - 1) {  // Prevent the last permutation from occurring.
            this->alpha().transformONVToNextPermutation(onv_alpha);
        }
    }  // Alpha address (Ia) loop.

    return diagonal;
}


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
VectorX<double> SpinResolvedONVBasis::evaluateOperatorMatrixVectorProduct(const ScalarRSQOneElectronOperator<double>& f, const VectorX<double>& x) const {

    if (f.numberOfOrbitals() != this->alpha().numberOfOrbitals()) {
        throw std::invalid_argument("SpinResolvedONVBasis::evaluateOperatorMatrixVectorProduct(const ScalarRSQOneElectronOperator<double>&, const VectorX<double>&): The number of orbitals of the spin-resolved ONV basis and the operator are incompatible.");
    }

    // Prepare some variables.
    const auto& alpha_couplings = this->alphaCouplings();

    const auto dim_alpha = static_cast<long>(this->alpha().dimension());  // Casting is required because of Eigen.
    const auto dim_beta = static_cast<long>(this->beta().dimension());


    // We can calculate the evaluation using re-mapped approach. We first map x as a dense matrix instead of a vector, and prepare a zero-initialized vector for storing the result.
    Eigen::Map<const Eigen::MatrixXd> x_map {x.data(), dim_beta, dim_alpha};
    VectorX<double> matvec = VectorX<double>::Zero(this->dimension());
    Eigen::Map<Eigen::MatrixXd> matvec_map {matvec.data(), dim_beta, dim_alpha};

    // The contributions can then be written very simply as matrix-matrix multiplications, taking advantage of mapped matvec representation. We use a sparse multiplication in order to reduce memory and speed impact.
    const auto H_a = this->alpha().evaluateOperatorSparse(f.alpha());
    const auto H_b = this->beta().evaluateOperatorSparse(f.beta());
    matvec_map += H_b * x_map + x_map * H_a;

    // We can safely return the vector representation of the matvec, because we have used Eigen's mapped representation to emplace its elements.
    return matvec;
}


/**
 *  Calculate the matrix-vector product of (the matrix representation of) a restricted two-electron operator with the given coefficient vector.
 *
 *  @param g                A restricted two-electron operator expressed in an orthonormal orbital basis.
 *  @param x                The coefficient vector of a linear expansion.
 *
 *  @return The coefficient vector of the linear expansion after being acted on with the given (matrix representation of) the two-electron operator.
 */
VectorX<double> SpinResolvedONVBasis::evaluateOperatorMatrixVectorProduct(const ScalarRSQTwoElectronOperator<double>& g, const VectorX<double>& x) const {

    // We choose to avoid code duplication by evaluating an equivalent restricted Hamiltonian with zero-valued core contributions.
    const auto zero = ScalarRSQOneElectronOperator<double>::Zero(g.numberOfOrbitals());
    const RSQHamiltonian<double> hamiltonian {zero, g};

    return this->evaluateOperatorMatrixVectorProduct(hamiltonian, x);
}


/**
 *  Calculate the matrix-vector product of (the matrix representation of) a restricted Hamiltonian with the given coefficient vector.
 *
 *  @param hamiltonian      A restricted Hamiltonian expressed in an orthonormal orbital basis.
 *  @param x                The coefficient vector of a linear expansion.
 *
 *  @return The coefficient vector of the linear expansion after being acted on with the given (matrix representation of) the Hamiltonian.
 */
VectorX<double> SpinResolvedONVBasis::evaluateOperatorMatrixVectorProduct(const RSQHamiltonian<double>& hamiltonian, const VectorX<double>& x) const {

    // We can avoid code duplication by delegating this method to the evaluation of an unrestricted Hamiltonian. This entails a small speed decrease, but doesn't change the order of the scaling.
    const auto h_unrestricted = ScalarUSQOneElectronOperator<double>::FromRestricted(hamiltonian.core());
    const auto g_unrestricted = ScalarUSQTwoElectronOperator<double>::FromRestricted(hamiltonian.twoElectron());
    const USQHamiltonian<double> unrestricted_hamiltonian {h_unrestricted, g_unrestricted};

    return this->evaluateOperatorMatrixVectorProduct(unrestricted_hamiltonian, x);
}


/**
 *  Calculate the matrix-vector product of (the matrix representation of) a Hubbard Hamiltonian with the given coefficient vector.
 *
 *  @param hamiltonian      A Hubbard Hamiltonian expressed in an orthonormal orbital basis.
 *  @param x                The coefficient vector of a linear expansion.
 *
 *  @return The coefficient vector of the linear expansion after being acted on with the given (matrix representation of) the Hamiltonian.
 */
VectorX<double> SpinResolvedONVBasis::evaluateOperatorMatrixVectorProduct(const HubbardHamiltonian<double>& hamiltonian, const VectorX<double>& x) const {

    if (hamiltonian.numberOfLatticeSites() != this->alpha().numberOfOrbitals()) {
        throw std::invalid_argument("SpinResolvedONVBasis::evaluateOperatorMatrixVectorProduct(const HubbardHamiltonian<double>&, const VectorX<double>&): The number of spatial orbitals of this ONV basis and the number of lattice sites for the Hubbard Hamiltonian are incompatible.");
    }


    // Prepare some variables.
    const auto dim_alpha = static_cast<long>(this->alpha().dimension());  // Casting is required because of Eigen.
    const auto dim_beta = static_cast<long>(this->beta().dimension());
    const auto diagonal = this->evaluateOperatorDiagonal(hamiltonian);
    const auto h = hamiltonian.core();


    // Calculate the Hubbard matrix-vector product, which is the sum of the alpha- and beta one-electron contributions plus the diagonal that contains the two-electron contributions.

    // We can calculate the evaluation using re-mapped approach. We first map x as a dense matrix instead of a vector, and calculate an initial value for the vector (from the diagonal two-electron evaluations) for storing the result.
    Eigen::Map<const Eigen::MatrixXd> x_map {x.data(), dim_beta, dim_alpha};
    VectorX<double> matvec = diagonal.cwiseProduct(x);
    Eigen::Map<Eigen::MatrixXd> matvec_map {matvec.data(), dim_beta, dim_alpha};

    // The contributions can then be written very simply as matrix-matrix multiplications, taking advantage of mapped matvec representation. We use a sparse multiplication in order to reduce memory and speed impact.
    const auto H_a = this->alpha().evaluateOperatorSparse(h.alpha());
    const auto H_b = this->alpha().evaluateOperatorSparse(h.beta());
    matvec_map += H_b * x_map + x_map * H_a;

    // We can safely return the vector representation of the matvec, because we have used Eigen's mapped representation to emplace its elements.
    return matvec;
}


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
VectorX<double> SpinResolvedONVBasis::evaluateOperatorDiagonal(const USQHamiltonian<double>& hamiltonian) const {

    if (hamiltonian.numberOfOrbitals() != this->alpha().numberOfOrbitals()) {
        throw std::invalid_argument("SpinResolvedONVBasis::evaluateOperatorDense(const USQHamiltonian<double>&): The number of orbitals of this ONV basis and the given Hamiltonian are incompatible.");
    }


    // Set up some variables.
    const auto dim_alpha = this->alpha().dimension();
    const auto dim_beta = this->beta().dimension();

    // We're going to use the effective one-electron contributions to simplify the algorithm.
    const auto k_a = (hamiltonian.core().alpha() + hamiltonian.twoElectron().alphaAlpha().effectiveOneElectronPartition()).parameters();
    const auto k_b = (hamiltonian.core().beta() + hamiltonian.twoElectron().betaBeta().effectiveOneElectronPartition()).parameters();

    const auto& g_a = hamiltonian.twoElectron().alphaAlpha().parameters();
    const auto& g_b = hamiltonian.twoElectron().betaBeta().parameters();
    auto const& g_ab = hamiltonian.twoElectron().alphaBeta().parameters();  // For contributions related to g_ab, we can use the relation g_ab(pqrs) = g_ba(rspq).


    VectorX<double> diagonal = VectorX<double>::Zero(this->dimension());


    SpinUnresolvedONV onv_alpha = this->alpha().constructONVFromAddress(0);
    SpinUnresolvedONV onv_beta = this->beta().constructONVFromAddress(0);
    for (size_t Ia = 0; Ia < dim_alpha; Ia++) {  // Ia loops over addresses of alpha ONVs

        this->beta().transformONVCorrespondingToAddress(onv_beta, 0);  // Reset the beta ONV.
        for (size_t Ib = 0; Ib < dim_beta; Ib++) {                     // Ib loops over addresses of beta ONVs

            for (size_t e_a = 0; e_a < this->alpha().numberOfElectrons(); e_a++) {  // loop over alpha electrons

                size_t p = onv_alpha.occupationIndexOf(e_a);
                diagonal(Ia * dim_beta + Ib) += k_a(p, p);

                for (size_t q = 0; q < this->alpha().numberOfOrbitals(); q++) {
                    if (onv_alpha.isOccupied(q)) {  // q is in Ia
                        diagonal(Ia * dim_beta + Ib) += 0.5 * g_a(p, p, q, q);
                    } else {  // q is not in I_alpha
                        diagonal(Ia * dim_beta + Ib) += 0.5 * g_a(p, q, q, p);
                    }

                    if (onv_beta.isOccupied(q)) {  // q is in Ib
                        diagonal(Ia * dim_beta + Ib) += g_ab(p, p, q, q);
                    }
                }  // q loop
            }      // e_a loop

            for (size_t e_b = 0; e_b < this->beta().numberOfElectrons(); e_b++) {  // loop over beta electrons

                size_t p = onv_beta.occupationIndexOf(e_b);
                diagonal(Ia * dim_beta + Ib) += k_b(p, p);

                for (size_t q = 0; q < this->alpha().numberOfOrbitals(); q++) {
                    if (onv_beta.isOccupied(q)) {  // q is in Ib
                        diagonal(Ia * dim_beta + Ib) += 0.5 * g_b(p, p, q, q);
                    } else {  // q is not in I_beta
                        diagonal(Ia * dim_beta + Ib) += 0.5 * g_b(p, q, q, p);
                    }
                }  // q loop
            }      // e_b loop

            if (Ib < dim_beta - 1) {  // prevent the last permutation from occurring
                this->beta().transformONVToNextPermutation(onv_beta);
            }
        }  // beta address (Ib) loop

        if (Ia < dim_alpha - 1) {  // prevent the last permutation from occurring
            this->alpha().transformONVToNextPermutation(onv_alpha);
        }
    }  // alpha address (Ia) loop

    return diagonal;
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
VectorX<double> SpinResolvedONVBasis::evaluateOperatorMatrixVectorProduct(const USQHamiltonian<double>& hamiltonian, const VectorX<double>& x) const {

    if (hamiltonian.numberOfOrbitals() != this->alpha().numberOfOrbitals()) {
        throw std::invalid_argument("SpinResolvedONVBasis::evaluateOperatorDense(const USQHamiltonian<double>&): The number of orbitals of this ONV basis and the given Hamiltonian are incompatible.");
    }

    // Prepare some variables.
    const auto& alpha_couplings = this->alphaCouplings();

    const auto dim_alpha = static_cast<long>(this->alpha().dimension());  // Casting is required because of Eigen.
    const auto dim_beta = static_cast<long>(this->beta().dimension());

    // In order to call the semantically correct APIs in the remainder of this method, we'll have to convert the pure alpha and pure beta part of the unrestricted Hamiltonian into a generalized representation.
    const auto& h_a = ScalarGSQOneElectronOperator<double>::FromUnrestrictedComponent(hamiltonian.core().alpha());
    const auto& g_aa = ScalarGSQTwoElectronOperator<double>::FromUnrestrictedComponent(hamiltonian.twoElectron().alphaAlpha());
    const GSQHamiltonian<double> alpha_hamiltonian {h_a, g_aa};

    const auto& h_b = ScalarGSQOneElectronOperator<double>::FromUnrestrictedComponent(hamiltonian.core().beta());
    const auto& g_bb = ScalarGSQTwoElectronOperator<double>::FromUnrestrictedComponent(hamiltonian.twoElectron().betaBeta());
    const GSQHamiltonian<double> beta_hamiltonian {h_b, g_bb};

    auto const& g_mixed = hamiltonian.twoElectron().alphaBeta();


    // We can calculate the 'pure spin evaluations', i.e. those only resulting exclusively from the alph and beta part, using re-mapped approach. We first map x as a dense matrix instead of a vector, and prepare a zero-initialized vector for storing the result.
    Eigen::Map<const Eigen::MatrixXd> x_map {x.data(), dim_beta, dim_alpha};
    VectorX<double> matvec = VectorX<double>::Zero(this->dimension());
    Eigen::Map<Eigen::MatrixXd> matvec_map {matvec.data(), dim_beta, dim_alpha};

    // The 'pure spin contributions' can then be written very simply as matrix-matrix multiplications, taking advantage of mapped matvec representation. We use a sparse multiplication in order to reduce memory and speed impact.
    const auto H_a = this->alpha().evaluateOperatorSparse(alpha_hamiltonian);
    const auto H_b = this->beta().evaluateOperatorSparse(beta_hamiltonian);

    matvec_map += H_b * x_map + x_map * H_a;


    // For the 'mixed spin contributions', i.e. those resulting from the alpha-beta (and beta-alpha) part of the two-electron part of the Hamiltonian, we can use the intermediate variables 'sigma' and 'theta'.
    for (size_t p = 0; p < this->alpha().numberOfOrbitals(); p++) {

        const auto& P = this->calculateOneElectronPartition(p, p, g_mixed);
        const auto& beta_two_electron_intermediate = this->beta().evaluateOperatorDense(P);

        // sigma(pp) * X * theta(pp)
        matvec_map += beta_two_electron_intermediate * (x_map * alpha_couplings[p * (this->alpha().numberOfOrbitals() + this->alpha().numberOfOrbitals() + 1 - p) / 2]);
        for (size_t q = p + 1; q < this->alpha().numberOfOrbitals(); q++) {

            const auto& P = this->calculateOneElectronPartition(p, q, g_mixed);
            const auto& beta_two_electron_intermediate = this->beta().evaluateOperatorDense(P);

            // (sigma(pq) + sigma(qp)) * X * theta(pq)
            matvec_map += beta_two_electron_intermediate * (x_map * alpha_couplings[p * (this->alpha().numberOfOrbitals() + this->alpha().numberOfOrbitals() + 1 - p) / 2 + q - p]);
        }
    }

    // We can safely return the vector representation of the matvec, because we have used Eigen's mapped representation to emplace its elements.
    return matvec;
}


}  // namespace GQCP
