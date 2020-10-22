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


// #include "ONVBasis/SpinResolvedONVBasis.hpp"
// #include "Operator/SecondQuantized/ModelHamiltonian/HubbardHamiltonian.hpp"


// namespace GQCP {


// /**
//  *  Hubbard builds a Hubbard Hamiltonian matrix in the full spin-resolved ONV basis
//  *
//  *  Hubbard distinguishes itself from FCI by explicitly implementing simplified matrix construction and matrix-vector implementations.
//  */
// class Hubbard {
// private:
//     SpinResolvedONVBasis onv_basis;  // the full spin-resolved ONV basis


// public:
//     // CONSTRUCTORS

//     /**
//      *  @param onv_basis       the full spin-resolved ONV basis
//      */
//     explicit Hubbard(const SpinResolvedONVBasis& onv_basis);


//     // DESTRUCTOR
//     ~Hubbard() = default;


//     // PUBLIC METHODS

//     /**
//      *  @param hubbard_hamiltonian              the Hubbard model Hamiltonian
//      *
//      *  @return the diagonal of the matrix representation of the Hubbard model Hamiltonian
//      */
//     VectorX<double> calculateDiagonal(const HubbardHamiltonian<double>& hubbard_hamiltonian) const;

//     /**
//      *  @param hubbard_hamiltonian              the Hubbard model Hamiltonian
//      *
//      *  @return the Hubbard Hamiltonian matrix
//      */
//     SquareMatrix<double> constructHamiltonian(const HubbardHamiltonian<double>& hubbard_hamiltonian) const;

//     /**
//      *  @param hubbard_hamiltonian              the Hubbard model Hamiltonian
//      *  @param x                                the vector upon which the Hubbard Hamiltonian acts
//      *  @param diagonal                         the diagonal of the Hubbard Hamiltonian matrix
//      *
//      *  @return the action of the Hubbard model Hamiltonian on the coefficient vector
//      */
//     VectorX<double> matrixVectorProduct(const HubbardHamiltonian<double>& hubbard_hamiltonian, const VectorX<double>& x, const VectorX<double>& diagonal) const;
// };


// }  // namespace GQCP
