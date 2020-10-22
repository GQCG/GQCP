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
 *  CONSTRUCTORS
 */

/**
 *  @param K            the number of spatial orbitals
 *  @param N_P          the number of electron pairs
 */
SeniorityZeroONVBasis::SeniorityZeroONVBasis(const size_t K, const size_t N_P) :
    K {K},
    N_P {N_P},
    dim {SeniorityZeroONVBasis::calculateDimension(K, N_P)} {}


/*
 *  STATIC PUBLIC METHODS
 */

/**
 *  @param K            the number of spatial orbitals
 *  @param N_P          the number of electron pairs
 * 
 *  @return the dimension of a seniority-zero ONV basis with the given number of spatial orbitals and electron pairs
 */
size_t SeniorityZeroONVBasis::calculateDimension(const size_t K, const size_t N_P) {

    return SpinUnresolvedONVBasis::calculateDimension(K, N_P);
}


/*
 *  PUBLIC METHODS
 */

/**
 *  Evaluate the diagonal of the operator
 *
 *  @param one_op               the one-electron operator in an orthonormal orbital basis to be evaluated in the ONV basis
 *
 *  @return the operator's diagonal evaluation in a vector with the dimension of the ONV basis
 */
// VectorX<double> SeniorityZeroONVBasis::evaluateOperatorDiagonal(const ScalarRSQOneElectronOperator<double>& one_op) const {

//     // Check if the argument is compatible.
//     const auto K = one_op.numberOfOrbitals();  // number of spatial orbitals

//     if (K != this->numberOfSpatialOrbitals()) {
//         throw std::invalid_argument("SeniorityZeroONVBasis::evaluateOperatorDiagonal(const ScalarRSQOneElectronOperator<double>&): The number of spatial orbitals for the ONV basis and one-electron operator are incompatible.");
//     }

//     // Prepare some variables to be used in the algorithm.
//     const auto dim = this->dimension();
//     const auto& f = one_op.parameters();

//     VectorX<double> diagonal = VectorX<double>::Zero(dim);


//     // Iterate over every proxy doubly-occupied ONV. Since we are actually using spin-unresolved ONVs, we should multiply contributions by 2.
//     this->forEach([&diagonal, &f](const SpinUnresolvedONV& onv, const size_t I) {
//         double value = 0;  // to be added to the diagonal

//         // Loop over every occupied orbital index and add the contribution.
//         onv.forEach([&value, &f](const size_t p) {
//             value += 2 * f(p, p);  // *2 because of seniority-zero
//         });

//         diagonal(I) += value;
//     });

//     return diagonal;
// }


// /**
//  *  Evaluate the diagonal of the matrix representation of a two-electron operator inside this seniority-zero ONV basis.
//  *
//  *  @param two_op               a two-electron operator expressed in an orthonormal orbital basis
//  *
//  *  @return the diagonal of the matrix representation of the two-electron operator in this seniority-zero ONV basis
//  */
// VectorX<double> SeniorityZeroONVBasis::evaluateOperatorDiagonal(const ScalarRSQTwoElectronOperator<double>& two_op) const {

//     // Check if the argument is compatible.
//     const auto K = two_op.numberOfOrbitals();  // number of spatial orbitals

//     if (K != this->numberOfSpatialOrbitals()) {
//         throw std::invalid_argument("SeniorityZeroONVBasis::evaluateOperatorDiagonal(const ScalarRSQOneElectronOperator<double>&): The number of spatial orbitals for the ONV basis and one-electron operator are incompatible.");
//     }

//     // Prepare some variables to be used in the algorithm.
//     const auto dim = this->dimension();
//     const auto& g = two_op.parameters();

//     VectorX<double> diagonal = VectorX<double>::Zero(dim);


//     // Iterate over every proxy doubly-occupied ONV. Since we are actually using spin-unresolved ONVs, we should multiply contributions by 2.
//     this->forEach([&diagonal, &g](const SpinUnresolvedONV& onv, const size_t I) {
//         double value = 0;  // to be added to the diagonal

//         // Loop over every occupied spinor index and add the contributions.
//         onv.forEach([&value, &g](const size_t p) {
//             value += g(p, p, p, p);  // 1/2*2 because of seniority-zero
//         });

//         // Loop over every pair of occupied spinor indices and add the contributions.
//         onv.forEach([&value, &g](const size_t p, const size_t q) {
//             // Since we are doing a restricted summation (p > q), we should multiply by 2 since the summand argument is symmetric upon interchanging p and q.
//             value += 2 * (2 * g(p, p, q, q) - g(p, q, q, p));
//         });

//         diagonal(I) += value;
//     });

//     return diagonal;
// }


/**
 *  Evaluate the diagonal of the matrix representation of a Hamiltonian inside this seniority-zero ONV basis.
 *
 *  @param sq_hamiltonian               a Hamiltonian expressed in an orthonormal orbital basis
 *
 *  @return the diagonal of the matrix representation of the Hamiltonian in this seniority-zero ONV basis
 */
// VectorX<double> SeniorityZeroONVBasis::evaluateOperatorDiagonal(const RSQHamiltonian<double>& sq_hamiltonian) const {

//     // We don't just use the sum of the one- and two-electron operator's diagonal representation because that would mean 2 iterations over the dimension of the ONV basis.


//     // Check if the argument is compatible.
//     const auto K = sq_hamiltonian.numberOfOrbitals();  // number of spatial orbitals

//     if (K != this->numberOfSpatialOrbitals()) {
//         throw std::invalid_argument("SeniorityZeroONVBasis::evaluateOperatorDiagonal(const RSQHamiltonian<double>&): The number of spatial orbitals for the ONV basis and one-electron operator are incompatible.");
//     }

//     // Prepare some variables to be used in the algorithm.
//     const auto dim = this->dimension();

//     const auto& h = sq_hamiltonian.core().parameters();
//     const auto& g = sq_hamiltonian.twoElectron().parameters();

//     VectorX<double> diagonal = VectorX<double>::Zero(dim);


//     // Iterate over every proxy doubly-occupied ONV. Since we are actually using spin-unresolved ONVs, we should multiply contributions by 2.
//     this->forEach([&diagonal, &h, &g](const SpinUnresolvedONV& onv, const size_t I) {
//         double value = 0;  // to be added to the diagonal

//         // Loop over every occupied spinor index and add the contributions.
//         onv.forEach([&value, &h, &g](const size_t p) {
//             value += 2 * h(p, p);    // *2 because of seniority zero
//             value += g(p, p, p, p);  // 1/2*2 because of seniority zero
//         });

//         // Loop over every pair of occupied spinor indices and add the contributions.
//         onv.forEach([&value, &g](const size_t p, const size_t q) {
//             // Since we are doing a restricted summation (p > q), we should multiply by 2 since the summand argument is symmetric upon interchanging p and q.
//             value += 2 * (2 * g(p, p, q, q) - g(p, q, q, p));
//         });

//         diagonal(I) += value;
//     });

//     return diagonal;
// }


/**
 *  Evaluate a one electron operator in a matrix vector product
 *
 *  @param one_op                       the one electron operator expressed in an orthonormal basis
 *  @param x                            the vector upon which the evaluation acts 
 *  @param diagonal                     the diagonal evaluated in the ONV basis
 *
 *  @return the one electron operator's matrix vector product in a vector with the dimensions of the ONV basis
 */
// VectorX<double> SeniorityZeroONVBasis::evaluateOperatorMatrixVectorProduct(const ScalarRSQOneElectronOperator<double>& one_op, const VectorX<double>& x, const VectorX<double>& diagonal) const {

//     const SpinResolvedSelectedONVBasis selected_onv_basis {*this};
//     return selected_onv_basis.evaluateOperatorMatrixVectorProduct(one_op, x, diagonal);
// }


/**
 *  Iterate over every (proxy) spin-resolved ONV in this seniority-zero ONV basis and apply the given callback.
 * 
 *  @param callback             a function to be called on every step during the iteration overall the ONVs. The arguments of the callback are the ONV and its address in this ONV basis.
 */
void SeniorityZeroONVBasis::forEach(const std::function<void(const SpinUnresolvedONV&, size_t)>& callback) const {

    // Create the first doubly-occupied ONV. Since in DOCI, alpha == beta, we can use the proxy ONV basis to treat them as one.
    const auto proxy_onv_basis = this->proxy();
    SpinUnresolvedONV onv = proxy_onv_basis.constructONVFromAddress(0);  // ONV with address 0

    for (size_t I = 0; I < dim; I++) {  // I loops over addresses of spin strings

        callback(onv, I);

        if (I < dim - 1) {  // prevent the last permutation from occurring
            proxy_onv_basis.transformONVToNextPermutation(onv);
        }
    }  // address (I) loop
}


}  // namespace GQCP
