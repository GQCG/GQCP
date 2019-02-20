// This file is part of GQCG-gqcp.
//
// Copyright (C) 2017-2019  the GQCG developers
//
// GQCG-gqcp is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// GQCG-gqcp is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-gqcp.  If not, see <http://www.gnu.org/licenses/>.
//
#include "HamiltonianBuilder/FrozenCore.hpp"

#include "utilities/linalg.hpp"
#include <utility>


namespace GQCP {

/*
 *  CONSTRUCTORS
 */

/**
 *  @param hamiltonian_builder           shared pointer to active (non-frozen core) Hamiltonian builder
 *  @param X                             the number of frozen orbitals
 */
FrozenCore::FrozenCore(std::shared_ptr<GQCP::HamiltonianBuilder> hamiltonian_builder, size_t X) :
        HamiltonianBuilder(),
        hamiltonian_builder (std::move(hamiltonian_builder)),
        X (X)
{}



/*
 *  OVERRIDDEN PUBLIC METHODS
 */

/**
 *  @param hamiltonian_parameters       the Hamiltonian parameters in an orthonormal orbital basis
 *
 *  @return the frozen core Hamiltonian matrix
 */
Eigen::MatrixXd FrozenCore::constructHamiltonian(const HamiltonianParameters& hamiltonian_parameters) const {

    // Freeze Hamiltonian parameters
    HamiltonianParameters frozen_ham_par = this->freezeHamiltonianParameters(hamiltonian_parameters, X);

    // calculate Hamiltonian matrix through conventional CI
    Eigen::MatrixXd total_hamiltonian = this->hamiltonian_builder->constructHamiltonian(frozen_ham_par);

    // diagonal correction
    Eigen::VectorXd diagonal = Eigen::VectorXd::Ones(this->get_fock_space()->get_dimension());
    auto frozen_core_diagonal = this->calculateFrozenCoreDiagonal(hamiltonian_parameters, this->X);
    total_hamiltonian += frozen_core_diagonal.asDiagonal();

    return total_hamiltonian;
}


/**
 *  @param hamiltonian_parameters       the Hamiltonian parameters in an orthonormal orbital basis
 *  @param x                            the vector upon which the Hamiltonian acts
 *  @param diagonal                     the diagonal of the Hamiltonian matrix
 *
 *  @return the action of the frozen core Hamiltonian on the coefficient vector
 */
Eigen::VectorXd FrozenCore::matrixVectorProduct(const HamiltonianParameters& hamiltonian_parameters, const Eigen::VectorXd& x, const Eigen::VectorXd& diagonal) const {

    HamiltonianParameters frozen_ham_par = this->freezeHamiltonianParameters(hamiltonian_parameters, X);

    // perform matvec in the active space with "frozen" Hamiltonian parameters
    return this->hamiltonian_builder->matrixVectorProduct(frozen_ham_par, x, diagonal);
}


/**
 *  @param hamiltonian_parameters       the Hamiltonian parameters in an orthonormal orbital basis
 *
 *  @return the diagonal of the matrix representation of the frozen core Hamiltonian
 */
Eigen::VectorXd FrozenCore::calculateDiagonal(const HamiltonianParameters& hamiltonian_parameters) const {

    HamiltonianParameters frozen_ham_par = this->freezeHamiltonianParameters(hamiltonian_parameters, this->X);

    // calculate diagonal in the active space with the "frozen" Hamiltonian parameters
    Eigen::VectorXd diagonal = this->hamiltonian_builder->calculateDiagonal(frozen_ham_par);

    // calculate diagonal for the frozen orbitals
    auto frozen_core_diagonal = this->calculateFrozenCoreDiagonal(hamiltonian_parameters, this->X);

    return diagonal + frozen_core_diagonal;
}



/*
 *  PUBLIC METHODS
 */

/**
 *  @param hamiltonian_parameters               the Hamiltonian parameters in an orthonormal orbital basis
 *  @param X                                    the number of frozen orbitals
 *
 *  @return a set of 'frozen' Hamiltonian parameters which cover two-electron integral evaluations from the active and inactive orbitals
 *  (see https://drive.google.com/file/d/1Fnhv2XyNO9Xw9YDoJOXU21_6_x2llntI/view?usp=sharing)
 */
HamiltonianParameters FrozenCore::freezeHamiltonianParameters(const HamiltonianParameters &hamiltonian_parameters, size_t X) const {

    size_t K = hamiltonian_parameters.get_K();
    size_t Kn = K - X;  // number of non-frozen orbitals

    std::shared_ptr<AOBasis> ao_basis;  // nullptr

    // copy one electron integrals from the non-frozen orbitals
    OneElectronOperator S (hamiltonian_parameters.get_S().get_matrix_representation().block(X, X, Kn, Kn));
    Eigen::MatrixXd h (hamiltonian_parameters.get_h().get_matrix_representation().block(X, X, Kn, Kn));

    const auto& g = hamiltonian_parameters.get_g();

    // copy two electron integrals from the non-frozen orbitals
    Eigen::Tensor<double, 4> g_SO = tensorBlockCreation(g.get_matrix_representation(), X, X, X, X);

    // Modify one electron integrals with the missing frozen orbital two electron integral evaluations according to the frozen core CI algorithm
    for (size_t i = 0; i < Kn; i++) {

        size_t q = i + X;  // map index to the non-modified Hamiltonian parameters

        for (size_t l = 0; l < X; l++) {  // iterate over the frozen orbitals
            h(i,i) += g(q, q, l, l);
            h(i,i) += g(l, l, q, q);
            h(i,i) -= g(q, l, l, q)/2;
            h(i,i) -= g(l, q, q, l)/2;
        }

        for (size_t j = i+1; j < Kn; j++) {

            size_t p = j + X;  // map index to the non-modified Hamiltonian parameters

            for (size_t l = 0; l < X; l++) {  // iterate over the frozen orbitals

                h(i,j) += g(q, p, l, l);
                h(i,j) += g(l, l, q, p);
                h(i,j) -= g(q, l, l, p)/2;
                h(i,j) -= g(l, p, q, l)/2;

                h(j,i) += g(p, q, l, l);
                h(j,i) += g(l, l, p, q);
                h(j,i) -= g(p, l, l, q)/2;
                h(j,i) -= g(l, q, p, l)/2;
            }
        }
    }

    OneElectronOperator H_new (h);
    Eigen::MatrixXd C = Eigen::MatrixXd(hamiltonian_parameters.get_T_total().block(X, X, Kn, Kn));
    TwoElectronOperator g_new (g_SO);

    return HamiltonianParameters(ao_basis, S, H_new, g_new, C);
}


/**
 *  @param hamiltonian_parameters              the Hamiltonian parameters in an orthonormal orbital basis
 *  @param X                                   the number of frozen orbitals
 *
 *  @return the diagonal from strictly evaluating the frozen orbitals in the Fock space
 */
Eigen::VectorXd FrozenCore::calculateFrozenCoreDiagonal(const HamiltonianParameters& hamiltonian_parameters, size_t X) const {

    const auto& g = hamiltonian_parameters.get_g();
    const auto& h = hamiltonian_parameters.get_h();

    double value = 0;

    // calculate the diagonal attributed to the frozen orbitals
    for (size_t i = 0; i < X; i++) {
        value += 2*h(i,i) + g(i,i,i,i);
        for (size_t j = i+1; j < X; j++) {
            value += 2*g(i,i,j,j);
            value += 2*g(j,j,i,i);
            value -= g(j,i,i,j);
            value -= g(i,j,j,i);
        }
    }

    Eigen::VectorXd diagonal = Eigen::VectorXd::Ones(this->get_fock_space()->get_dimension());
    return value * diagonal;
}


}  // namespace GQCP

