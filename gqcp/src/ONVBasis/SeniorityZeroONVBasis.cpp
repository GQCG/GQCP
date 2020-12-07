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

#include "ONVBasis/SeniorityZeroONVBasis.hpp"

#include "ONVBasis/SpinResolvedSelectedONVBasis.hpp"
#include "ONVBasis/SpinUnresolvedONVBasis.hpp"


namespace GQCP {


/*
 *  MARK: Constructors
 */

/**
 *  @param K            The number of spatial orbitals.
 *  @param N_P          The number of electron pairs.
 */
SeniorityZeroONVBasis::SeniorityZeroONVBasis(const size_t K, const size_t N_P) :
    K {K},
    N_P {N_P} {}


/*
 *  MARK: General information
 */


/**
 *  Calculate the dimension of a seniority-zero ONV basis with a given number of spatial orbitals and number of electron pairs.
 * 
 *  @param K            The number of spatial orbitals.
 *  @param N_P          The number of electron pairs.
 * 
 *  @return The dimension of a seniority-zero ONV basis.
 */
size_t SeniorityZeroONVBasis::calculateDimension(const size_t K, const size_t N_P) {

    return SpinUnresolvedONVBasis::calculateDimension(K, N_P);
}


/*
 *  MARK: Iterations
 */

/**
 *  Iterate over every (proxy) spin-(un)resolved ONV in this seniority-zero ONV basis and apply the given callback.
 * 
 *  @param callback             The function to be called on every step during the iteration overall the ONVs. The arguments of the callback are the ONV and its address in this ONV basis.
 */
void SeniorityZeroONVBasis::forEach(const std::function<void(const SpinUnresolvedONV&, size_t)>& callback) const {

    // Create the first doubly-occupied ONV. Since in a seniority-zero ONV basis, alpha == beta, we can use the proxy ONV basis to treat them as equal.
    const auto dim = this->dimension();
    const auto proxy_onv_basis = this->proxy();
    SpinUnresolvedONV onv = proxy_onv_basis.constructONVFromAddress(0);  // Initialize the ONV with address 0.
    for (size_t I = 0; I < dim; I++) {                                   // I loops over addresses of spin strings

        callback(onv, I);

        if (I < dim - 1) {  //Pprevent the last permutation from occurring.
            proxy_onv_basis.transformONVToNextPermutation(onv);
        }
    }  // Address (I) loop.
}


/*
 *  MARK: Dense restricted operator evaluations
 */

/**
 *  Calculate the dense matrix representation of a Hubbard Hamiltonian in this ONV basis.
 *
 *  @param hamiltonian      A Hubbard Hamiltonian expressed in an orthonormal orbital basis.
 *
 *  @return A dense matrix represention of the Hamiltonian.
 */
SquareMatrix<double> SeniorityZeroONVBasis::evaluateOperatorDense(const RSQHamiltonian<double>& hamiltonian) const {

    if (hamiltonian.numberOfOrbitals() != this->numberOfSpatialOrbitals()) {
        throw std::invalid_argument("SeniorityZeroONVBasis::evaluateOperatorDense(const RSQHamiltonian<double>&): The number of spatial orbitals for the ONV basis and Hamiltonian are incompatible.");
    }


    // Prepare some variables to be used in the algorithm.
    const size_t N_P = this->numberOfElectronPairs();
    const size_t dim = this->dimension();

    SquareMatrix<double> H = SquareMatrix<double>::Zero(dim);  // The matrix representation of the Hamiltonian.
    const auto diagonal = this->evaluateOperatorDiagonal(hamiltonian);
    const auto& g = hamiltonian.twoElectron().parameters();


    // Use a proxy ONV basis to treat alpha- and beta- ONVs as equal and multiply all contributions by 2.
    const auto proxy_onv_basis = this->proxy();
    auto onv = proxy_onv_basis.constructONVFromAddress(0);  // Create the ONV with address 0.
    for (size_t I = 0; I < dim; I++) {                      // I loops over all the addresses of the ONVs.

        H(I, I) += diagonal(I);

        for (size_t e1 = 0; e1 < N_P; e1++) {            // E1 (electron 1) loops over the number of electrons.
            const size_t p = onv.occupationIndexOf(e1);  // Retrieve the index of the orbital that the electron occupies.

            // Remove the weight from the initial address I, because we annihilate.
            size_t address = I - proxy_onv_basis.vertexWeight(p, e1 + 1);

            // The e2 iteration counts the number of encountered electrons for the creation operator.
            // We only consider greater addresses than the initial one (because of symmetry), hence we only count electron after the annihilated electron (e1).
            size_t e2 = e1 + 1;
            size_t q = p + 1;

            // perform a shift  TODO: clarify what shift.
            proxy_onv_basis.shiftUntilNextUnoccupiedOrbital<1>(onv, address, q, e2);

            while (q < K) {
                const size_t J = address + proxy_onv_basis.vertexWeight(q, e2);

                H(I, J) += g(p, q, p, q);
                H(J, I) += g(p, q, p, q);

                q++;  // Go to the next orbital.

                // Perform a shift;  TODO: clarify what shift.
                proxy_onv_basis.shiftUntilNextUnoccupiedOrbital<1>(onv, address, q, e2);
            }  // Creation.
        }      // E1 loop (annihilation).

        if (I < dim - 1) {  // Prevent the last permutation from occurring.
            proxy_onv_basis.transformONVToNextPermutation(onv);
        }
    }  // Address (I) loop.

    return H;
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
VectorX<double> SeniorityZeroONVBasis::evaluateOperatorDiagonal(const ScalarRSQOneElectronOperator<double>& f_op) const {

    if (f_op.numberOfOrbitals() != this->numberOfSpatialOrbitals()) {
        throw std::invalid_argument("SeniorityZeroONVBasis::evaluateOperatorDiagonal(const ScalarRSQOneElectronOperator<double>&): The number of spatial orbitals for the ONV basis and one-electron operator are incompatible.");
    }


    // Prepare some variables to be used in the algorithm.
    const auto dim = this->dimension();
    const auto& f = f_op.parameters();

    VectorX<double> diagonal = VectorX<double>::Zero(dim);


    // Iterate over every proxy doubly-occupied ONV. Since we are actually using spin-unresolved ONVs, we should multiply contributions by 2.
    this->forEach([&diagonal, &f](const SpinUnresolvedONV& onv, const size_t I) {
        double value = 0;  // to be added to the diagonal

        // Loop over every occupied orbital index and add the contribution.
        onv.forEach([&value, &f](const size_t p) {
            value += 2 * f(p, p);  //Factor  *2 because of seniority-zero.
        });

        diagonal(I) += value;
    });

    return diagonal;
}


/**
 *  Calculate the diagonal of the matrix representation of a restricted two-electron operator in this ONV basis.
 *
 *  @param g                A restricted two-electron operator expressed in an orthonormal orbital basis.
 *
 *  @return The diagonal of the dense matrix represention of the two-electron operator.
 */
VectorX<double> SeniorityZeroONVBasis::evaluateOperatorDiagonal(const ScalarRSQTwoElectronOperator<double>& g_op) const {

    // Check if the argument is compatible.
    if (g_op.numberOfOrbitals() != this->numberOfSpatialOrbitals()) {
        throw std::invalid_argument("SeniorityZeroONVBasis::evaluateOperatorDiagonal(const ScalarRSQOneElectronOperator<double>&): The number of spatial orbitals for the ONV basis and one-electron operator are incompatible.");
    }

    // Prepare some variables to be used in the algorithm.
    const auto dim = this->dimension();
    const auto& g = g_op.parameters();

    VectorX<double> diagonal = VectorX<double>::Zero(dim);


    // Iterate over every proxy doubly-occupied ONV. Since we are actually using spin-unresolved ONVs, we should multiply contributions by 2.
    this->forEach([&diagonal, &g](const SpinUnresolvedONV& onv, const size_t I) {
        double value = 0;  // to be added to the diagonal

        // Loop over every occupied spinor index and add the contributions.
        onv.forEach([&value, &g](const size_t p) {
            value += g(p, p, p, p);  // Factor 1/2*2 because of seniority-zero.
        });

        // Loop over every pair of occupied spinor indices and add the contributions.
        onv.forEach([&value, &g](const size_t p, const size_t q) {
            // Since we are doing a restricted summation (p > q), we should multiply by 2 since the summand argument is symmetric upon interchanging p and q.
            value += 2 * (2 * g(p, p, q, q) - g(p, q, q, p));
        });

        diagonal(I) += value;
    });

    return diagonal;
}


/**
 *  Calculate the diagonal of the dense matrix representation of a restricted Hamiltonian in this ONV basis.
 *
 *  @param hamiltonian      A restricted Hamiltonian expressed in an orthonormal orbital basis.
 *
 *  @return The diagonal of the dense matrix represention of the Hamiltonian.
 */
VectorX<double> SeniorityZeroONVBasis::evaluateOperatorDiagonal(const RSQHamiltonian<double>& hamiltonian) const {

    // We don't just use the sum of the one- and two-electron operator's diagonal representation because that would mean 2 iterations over the dimension of the ONV basis.


    // Check if the argument is compatible.
    const auto K = hamiltonian.numberOfOrbitals();  // The number of spatial orbitals.

    if (K != this->numberOfSpatialOrbitals()) {
        throw std::invalid_argument("SeniorityZeroONVBasis::evaluateOperatorDiagonal(const RSQHamiltonian<double>&): The number of spatial orbitals for the ONV basis and one-electron operator are incompatible.");
    }

    // Prepare some variables to be used in the algorithm.
    const auto dim = this->dimension();

    const auto& h = hamiltonian.core().parameters();
    const auto& g = hamiltonian.twoElectron().parameters();

    VectorX<double> diagonal = VectorX<double>::Zero(dim);


    // Iterate over every proxy doubly-occupied ONV. Since we are actually using spin-unresolved ONVs, we should multiply contributions by 2.
    this->forEach([&diagonal, &h, &g](const SpinUnresolvedONV& onv, const size_t I) {
        double value = 0;  // to be added to the diagonal

        // Loop over every occupied spinor index and add the contributions.
        onv.forEach([&value, &h, &g](const size_t p) {
            value += 2 * h(p, p);    // Factor *2 because of seniority zero.
            value += g(p, p, p, p);  // Factor 1/2*2 because of seniority zero.
        });

        // Loop over every pair of occupied spinor indices and add the contributions.
        onv.forEach([&value, &g](const size_t p, const size_t q) {
            // Since we are doing a restricted summation (p > q), we should multiply by 2 since the summand argument is symmetric upon interchanging p and q.
            value += 2 * (2 * g(p, p, q, q) - g(p, q, q, p));
        });

        diagonal(I) += value;
    });

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
VectorX<double> SeniorityZeroONVBasis::evaluateOperatorMatrixVectorProduct(const ScalarRSQOneElectronOperator<double>& f, const VectorX<double>& x) const {

    const SpinResolvedSelectedONVBasis selected_onv_basis {*this};
    return selected_onv_basis.evaluateOperatorMatrixVectorProduct(f, x);
}


/**
 *  Calculate the matrix-vector product of (the matrix representation of) a restricted Hamiltonian with the given coefficient vector.
 *
 *  @param hamiltonian      A restricted Hamiltonian expressed in an orthonormal orbital basis.
 *  @param x                The coefficient vector of a linear expansion.
 *
 *  @return The coefficient vector of the linear expansion after being acted on with the given (matrix representation of) the Hamiltonian.
 */
VectorX<double> SeniorityZeroONVBasis::evaluateOperatorMatrixVectorProduct(const RSQHamiltonian<double>& hamiltonian, const VectorX<double>& x) const {

    if (hamiltonian.numberOfOrbitals() != this->numberOfSpatialOrbitals()) {
        throw std::invalid_argument("DOCI::matrixVectorProduct(const RSQHamiltonian<double>&, const VectorX<double>&, const VectorX<double>&): The number of spatial orbitals for the ONV basis and Hamiltonian are incompatible.");
    }

    // Prepare some variables to be used in the algorithm.
    const size_t N_P = this->numberOfElectronPairs();
    const size_t dim = this->dimension();

    const auto& g = hamiltonian.twoElectron().parameters();


    // Initialize the resulting matrix-vector product from the diagonal contributions.
    VectorX<double> matvec = this->evaluateOperatorDiagonal(hamiltonian).cwiseProduct(x);

    // Create the first doubly-occupied ONV basis. Since in DOCI, alpha == beta, we can use the proxy ONV basis to treat them as one and multiply all contributions by 2.
    const auto proxy_onv_basis = this->proxy();
    auto onv = proxy_onv_basis.constructONVFromAddress(0);  // ONV with address 0.
    for (size_t I = 0; I < dim; I++) {                      // I loops over all the addresses of the ONV.

        // Using container values of type double reduce the number of times a vector has to be read from/written to.
        double value = 0;
        const double x_I = x(I);

        for (size_t e1 = 0; e1 < N_P; e1++) {            // E1 (electron 1) loops over the (number of) electrons.
            const size_t p = onv.occupationIndexOf(e1);  // Retrieve the index of a given electron.

            // Remove the weight from the initial address I, because we annihilate.
            size_t address = I - proxy_onv_basis.vertexWeight(p, e1 + 1);

            // The e2 iteration counts the number of encountered electrons for the creation operator.
            // We only consider greater addresses than the initial one (because of symmetry), hence we only count electron after the annihilated electron (e1).
            size_t e2 = e1 + 1;
            size_t q = p + 1;

            // Perform a shift.  TODO: clarify what shift.
            proxy_onv_basis.shiftUntilNextUnoccupiedOrbital<1>(onv, address, q, e2);

            while (q < K) {
                const size_t J = address + proxy_onv_basis.vertexWeight(q, e2);

                value += g(p, q, p, q) * x(J);
                matvec(J) += g(p, q, p, q) * x_I;

                q++;  // Go to the next orbital.

                // Perform a shift.  TODO: clarify what shift.
                proxy_onv_basis.shiftUntilNextUnoccupiedOrbital<1>(onv, address, q, e2);
            }  // Creation.
        }      // E1 loop (annihilation).

        if (I < dim - 1) {  // Prevent the last permutation.
            proxy_onv_basis.transformONVToNextPermutation(onv);
        }

        matvec(I) += value;
    }  // Address (I) loop.

    return matvec;
}


}  // namespace GQCP
