// // This file is part of GQCG-GQCP.
// //
// // Copyright (C) 2017-2020  the GQCG developers
// //
// // GQCG-GQCP is free software: you can redistribute it and/or modify
// // it under the terms of the GNU Lesser General Public License as published by
// // the Free Software Foundation, either version 3 of the License, or
// // (at your option) any later version.
// //
// // GQCG-GQCP is distributed in the hope that it will be useful,
// // but WITHOUT ANY WARRANTY; without even the implied warranty of
// // MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// // GNU Lesser General Public License for more details.
// //
// // You should have received a copy of the GNU Lesser General Public License
// // along with GQCG-GQCP.  If not, see <http://www.gnu.org/licenses/>.

// #pragma once


// #include "ONVBasis/SpinResolvedSelectedONVBasis.hpp"
// #include "QCMethod/CI/HamiltonianBuilder/HamiltonianBuilder.hpp"


// namespace GQCP {


// /**
//  *  Typedef for a type of function that handles where to 'put' a calculated value (e.g. matrix or vector).
//  *
//  *  @param I        index or address of an SpinUnresolvedONV
//  *  @param J        index or address of SpinUnresolvedONV that couples with SpinUnresolvedONV I
//  *  @param value    value related to the coupling
//  */
// using PassToMethod = std::function<void(size_t I, size_t J, double value)>;


// /**
//  *  SelectedCI builds a Hamiltonian matrix in a spin-resolved selected ONV basis
//  */
// class SelectedCI: public HamiltonianBuilder {
// private:
//     SpinResolvedSelectedONVBasis onv_basis;  // the spin-resolved selected ONV basis


// public:
//     // CONSTRUCTORS

//     /**
//      *  @param onv_basis               the spin-resolved selected ONV basis
//      */
//     explicit SelectedCI(const SpinResolvedSelectedONVBasis& onv_basis);


//     // DESTRUCTOR

//     /**
//      *  The default destructor.
//      */
//     ~SelectedCI() = default;


//     // PUBLIC OVERRIDDEN METHODS

//     /**
//      *  @param sq_hamiltonian               the Hamiltonian expressed in an orthonormal basis
//      *
//      *  @return the diagonal of the matrix representation of the SelectedCI Hamiltonian
//      */
//     VectorX<double> calculateDiagonal(const RSQHamiltonian<double>& sq_hamiltonian) const override;

//     /**
//      *  @param sq_hamiltonian               the Hamiltonian expressed in an orthonormal basis
//      *
//      *  @return the SelectedCI Hamiltonian matrix
//      */
//     SquareMatrix<double> constructHamiltonian(const RSQHamiltonian<double>& sq_hamiltonian) const override;

//     /**
//      *  @param sq_hamiltonian               the Hamiltonian expressed in an orthonormal basis
//      *  @param x                            the vector upon which the SelectedCI Hamiltonian acts
//      *  @param diagonal                     the diagonal of the SelectedCI Hamiltonian matrix
//      *
//      *  @return the action of the SelectedCI Hamiltonian on the coefficient vector
//      */
//     VectorX<double> matrixVectorProduct(const RSQHamiltonian<double>& sq_hamiltonian, const VectorX<double>& x, const VectorX<double>& diagonal) const override;

//     /**
//      *  @return the ONV basis that is associated with this HamiltonianBuilder
//      */
//     const BaseONVBasis* onvBasis() const override { return &onv_basis; }


//     // PUBLIC METHODS

//     /**
//      *  Evaluate all Hamiltonian elements, putting the results in the Hamiltonian matrix or matvec through the `method` function
//      *  This function is used both in `constructHamiltonian()` and `matrixVectorProduct()` to avoid duplicate code.
//      *
//      *  @param sq_hamiltonian           the Hamiltonian expressed in an orthonormal basis
//      *  @param method                   the method depending to how you wish to construct the Hamiltonian
//      */
//     void evaluateHamiltonianElements(const RSQHamiltonian<double>& sq_hamiltonian, const PassToMethod& method) const;
// };


// }  // namespace GQCP
