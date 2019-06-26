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
#include "Properties/AP1roGElectricalResponseSolver.hpp"

#include "Geminals/AP1roGPSEs.hpp"


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param G            the converged geminal coefficients
 *  @param lambda       the corresponding Lagrange multipliers
 */
AP1roGElectricalResponseSolver::AP1roGElectricalResponseSolver(const AP1roGGeminalCoefficients& G, const AP1roGVariables& lambda) :
    N_P (G.get_N_P()),
    G (G),
    lambda (lambda)
{}



/*
 *  PUBLIC OVERRIDDEN METHODS
 */

/**
 *  @param ham_par                  the Hamiltonian parameters
 * 
 *  @return the parameter response constant (k_p), i.e. the first-order parameter partial derivative of the PSEs function (i.e. the Jacobian of the PSEs)
 */
SquareMatrix<double> AP1roGElectricalResponseSolver::calculateParameterResponseConstant(const HamiltonianParameters<double>& ham_par) const {

    AP1roGPSEs pses (ham_par, this->N_P);
    const auto J = pses.evaluateJacobian(this->G);
    return SquareMatrix<double>(J.asMatrix());
}


/**
 *  @param dipole_integrals         the dipole integrals in an orthonormal orbital basis
 * 
 *  @return the parameter response force (F_p), i.e. the first-order perturbation derivative of the PSEs
 */
Matrix<double, Dynamic, 3> AP1roGElectricalResponseSolver::calculateParameterResponseForce(const std::array<OneElectronOperator<double>, 3>& dipole_integrals) const {

    const auto K = dipole_integrals[0].get_K();
    const auto dim = this->G.numberOfGeminalCoefficients(this->N_P, K);


    // Calculate every component separately
    Matrix<double, Dynamic, 3> F_p = Matrix<double, Dynamic, 3>::Zero(dim, 3);
    for (size_t m = 0; m < 3; m++) {
        
        const auto mu_m = dipole_integrals[m];

        // Calculate the m-th component of the parameter response force F_p
        BlockMatrix<double> F_p_m (0, this->N_P, this->N_P, K);
        for (size_t i = 0; i < this->N_P; i++) {
            for (size_t a = this->N_P; a < K; a++) {
                F_p_m(i,a) = 2 * (mu_m(i,i) - mu_m(a,a)) * this->G(i,a);
            }
        }


        // Place the component in the larger matrix
        F_p.col(m) = F_p_m.asVector();
    }


    return F_p;
}



/*
 *  PUBLIC METHODS
 */

/**
 *  @param ham_par                  the Hamiltonian parameters
 * 
 *  @return the Lagrangian multiplier response constant (k_lambda), i.e. the transpose of the parameter multiplier response constant
 */
SquareMatrix<double> AP1roGElectricalResponseSolver::calculateMultiplierResponseConstant(const HamiltonianParameters<double>& ham_par) const {

    return this->calculateParameterResponseConstant(ham_par).transpose();
}


/**
 *  @param dipole_integrals         the dipole integrals in an orthonormal orbital basis
 * 
 *  @return the explicit (i.e. the first part of the) Lagrangian multiplier response force, A_lambda
 */
Matrix<double, Dynamic, 3> AP1roGElectricalResponseSolver::calculateExplicitMultiplierResponseForce(const std::array<OneElectronOperator<double>, 3>& dipole_integrals) const {

    const auto K = dipole_integrals[0].get_K();
    const auto dim = this->G.numberOfGeminalCoefficients(this->N_P, K);

    // Calculate A_lambda by constructing its separate components
    Matrix<double, Dynamic, 3> A_lambda = Matrix<double, Dynamic, 3>::Zero(dim, 3);
    for (size_t m = 0; m < 3; m++) {

        const auto mu_m = dipole_integrals[m];

        // Calculate the m-th component
        BlockMatrix<double> A_lambda_m (0, this->N_P, this->N_P, K);
        for (size_t i = 0; i < this->N_P; i++) {
            for (size_t a = this->N_P; a < K; a++) {
                A_lambda_m(i,a) = 2 * this->lambda(i,a) * (mu_m(i,i) - mu_m(a,a));
            }
        }


        // Place the component in the larger matrix
        A_lambda.col(m) = A_lambda_m.asVector();
    }

    return A_lambda;
}


/**
 *  @param ham_par                  the Hamiltonian parameters
 *
 *  @return the multiplier force constant of the implicit part (i.e. the second part of the) Lagrangian multiplier response, B_lambda
 */
BlockRankFourTensor<double> AP1roGElectricalResponseSolver::calculateImplicitMultiplierResponseForceConstant(const HamiltonianParameters<double>& ham_par) const {

    const auto K = ham_par.get_K();
    const auto& g = ham_par.get_g();

    BlockRankFourTensor<double> B_lambda (0, this->N_P, this->N_P, K, 
                                          0, this->N_P, this->N_P, K);

    for (size_t i = 0; i < this->N_P; i++) {
        for (size_t a = this->N_P; a < K; a++) {
            for (size_t j = 0; j < this->N_P; j++) {
                for (size_t b = this->N_P; b < K; b++) {
                    double value {0.0};

                    value += this->lambda(i,b) * g(j,a,j,a) + this->lambda(j,a) * g(i,b,i,b);

                    if (i == j) {
                        value -= 2 * (this->lambda(i,b) * g(i,a,i,a) + this->lambda(i,a) * g(i,b,i,b));
                    }

                    if (a == b) {
                        value -= 2*(this->lambda(j,a) * g(i,a,i,a) + this->lambda(i,a) * g(j,a,j,a));
                    }

                    if ((i == j) && (a == b)) {
                        value += 4 * this->lambda(i,a) * g(i,a,i,a);
                    }
                }
            }
        }
    }

    return B_lambda;
}


/**
 *  @param ham_par                  the Hamiltonian parameters
 *  @param dipole_integrals         the dipole integrals in an orthonormal orbital basis
 *  @param x                        the first-order parameter response
 * 
 *  @return the Lagrangian multiplier response force (F_lambda)
 */
Matrix<double, Dynamic, 3> AP1roGElectricalResponseSolver::calculateMultiplierResponseForce(const HamiltonianParameters<double>& ham_par, const std::array<OneElectronOperator<double>, 3>& dipole_integrals, const Matrix<double, Dynamic, 3>& x) const {

    const auto A_lambda = this->calculateExplicitMultiplierResponseForce(dipole_integrals);
    const auto B_lambda = this->calculateImplicitMultiplierResponseForceConstant(ham_par);

    // Calculate F_lambda as F_lambda = A_lambda + B_lambda x
    return A_lambda + B_lambda.asMatrix() * x;
}


/**
 *  Solve the linear response equations for the Lagrangian multiplier response
 * 
 *  @param ham_par                  the Hamiltonian parameters
 *  @param dipole_integrals         the dipole integrals in an orthonormal orbital basis
 *  @param x                        the first-order parameter response
 * 
 *  @return the Lagrangian multiplier reponse y
 */
Matrix<double, Dynamic, 3> AP1roGElectricalResponseSolver::calculateMultiplierResponse(const HamiltonianParameters<double>& ham_par, const std::array<OneElectronOperator<double>, 3>& dipole_integrals, const Matrix<double, Dynamic, 3>& x) const {

    const auto k_lambda = this->calculateMultiplierResponseConstant(ham_par);
    const auto F_lambda = this->calculateMultiplierResponseForce(ham_par, dipole_integrals, x);


    // Solve k_lambda y = -F_lambda
    Eigen::HouseholderQR<Eigen::MatrixXd> linear_solver (k_lambda);
    const Matrix<double, Dynamic, 3> y = linear_solver.solve(-F_lambda);

    return y;
}


}  // namespace GQCP
