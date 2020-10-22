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

#include "Processing/Properties/vAP1roGElectricalResponseSolver.hpp"

#include "Basis/SpinorBasis/OrbitalSpace.hpp"
#include "Mathematical/Optimization/LinearEquation/LinearEquationEnvironment.hpp"
#include "Mathematical/Optimization/LinearEquation/LinearEquationSolver.hpp"


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param vap1rog          the optimal vAP1roG parameters
 */
vAP1roGElectricalResponseSolver::vAP1roGElectricalResponseSolver(const QCModel::vAP1roG& vap1rog) :
    vap1rog {vap1rog} {}


/*
 *  PUBLIC OVERRIDDEN METHODS
 */

/**
 *  @param sq_hamiltonian                   the Hamiltonian parameters
 * 
 *  @return the parameter response constant (k_p), i.e. the first-order parameter partial derivative of the PSEs, which is the Jacobian of the PSEs
 */
SquareMatrix<double> vAP1roGElectricalResponseSolver::calculateParameterResponseConstant(const RSQHamiltonian<double>& sq_hamiltonian) const {

    return QCModel::AP1roG(this->vap1rog.geminalCoefficients()).calculatePSEJacobian(sq_hamiltonian).asMatrix();
}


/**
 *  @param dipole_op                        the dipole integrals expressed in an orthonormal orbital basis
 * 
 *  @return the parameter response force (F_p) as an (Nx3)-matrix, i.e. the first-order perturbation derivative of the PSEs
 */
Matrix<double, Dynamic, 3> vAP1roGElectricalResponseSolver::calculateParameterResponseForce(const VectorRSQOneElectronOperator<double>& dipole_op) const {

    // Prepare some variables.
    const auto& G = this->vap1rog.geminalCoefficients();
    const auto dim = G.count();

    // Create an occupied-virtual orbital space.
    const auto K = G.numberOfSpatialOrbitals();
    const auto N_P = G.numberOfElectronPairs();
    const auto orbital_space = OrbitalSpace::Implicit({{OccupationType::k_occupied, N_P}, {OccupationType::k_virtual, K - N_P}});  // N_P occupied (spatial) orbitals, K-N_P virtual (spatial) orbitals


    // Calculate every component separately.
    Matrix<double, Dynamic, 3> F_p = Matrix<double, Dynamic, 3>::Zero(dim, 3);
    for (size_t m = 0; m < 3; m++) {

        const auto mu_m = dipole_op(m).parameters();

        // Calculate the m-th component of the parameter response force F_p.
        auto F_p_m = orbital_space.initializeRepresentableObjectFor<double>(OccupationType::k_occupied, OccupationType::k_virtual);  // create a representation for occupied-virtual objects
        for (const auto& i : orbital_space.indices(OccupationType::k_occupied)) {
            for (const auto& a : orbital_space.indices(OccupationType::k_virtual)) {
                F_p_m(i, a) = 2 * (mu_m(i, i) - mu_m(a, a)) * G(i, a);
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
 *  @param sq_hamiltonian                   the Hamiltonian expressed in an orthonormal orbital basis
 * 
 *  @return the Lagrangian multiplier response constant (k_lambda), which is the transpose of the parameter multiplier response constant
 */
SquareMatrix<double> vAP1roGElectricalResponseSolver::calculateMultiplierResponseConstant(const RSQHamiltonian<double>& sq_hamiltonian) const {

    return this->calculateParameterResponseConstant(sq_hamiltonian).transpose();
}


/**
 *  @param dipole_op                        the dipole integrals expressed in an orthonormal orbital basis
 * 
 *  @return the explicit (i.e. the first part of the) Lagrangian multiplier response force, A_lambda
 */
Matrix<double, Dynamic, 3> vAP1roGElectricalResponseSolver::calculateExplicitMultiplierResponseForce(const VectorRSQOneElectronOperator<double> dipole_op) const {

    // Prepare some variables.
    const auto& G = this->vap1rog.geminalCoefficients();
    const auto& lambda = this->vap1rog.lagrangeMultipliers();
    const auto dim = G.count();

    // Create an occupied-virtual orbital space.
    const auto K = G.numberOfSpatialOrbitals();
    const auto N_P = G.numberOfElectronPairs();
    const auto orbital_space = OrbitalSpace::Implicit({{OccupationType::k_occupied, N_P}, {OccupationType::k_virtual, K - N_P}});  // N_P occupied (spatial) orbitals, K-N_P virtual (spatial) orbitals


    // Calculate A_lambda by constructing its separate components.
    Matrix<double, Dynamic, 3> A_lambda = Matrix<double, Dynamic, 3>::Zero(dim, 3);
    for (size_t m = 0; m < 3; m++) {

        const auto mu_m = dipole_op(m).parameters();

        // Calculate the m-th component.
        auto A_lambda_m = orbital_space.initializeRepresentableObjectFor<double>(OccupationType::k_occupied, OccupationType::k_virtual);  // create a representation for occupied-virtual objects
        for (const auto& i : orbital_space.indices(OccupationType::k_occupied)) {
            for (const auto& a : orbital_space.indices(OccupationType::k_virtual)) {
                A_lambda_m(i, a) = 2 * lambda(i, a) * (mu_m(i, i) - mu_m(a, a));
            }
        }


        // Place the component in the larger matrix.
        A_lambda.col(m) = A_lambda_m.asVector();
    }

    return A_lambda;
}


/**
 *  @param sq_hamiltonian                   the Hamiltonian expressed in an orthonormal orbital basis
 *
 *  @return the multiplier force constant of the implicit part (i.e. the second part of the) Lagrangian multiplier response, B_lambda
 */
ImplicitRankFourTensorSlice<double> vAP1roGElectricalResponseSolver::calculateImplicitMultiplierResponseForceConstant(const RSQHamiltonian<double>& sq_hamiltonian) const {

    // Prepare some variables.
    const auto& G = this->vap1rog.geminalCoefficients();
    const auto& lambda = this->vap1rog.lagrangeMultipliers();
    const auto dim = G.count();

    // Prepare an occupied-virtual orbital space.
    const auto N_P = G.numberOfElectronPairs();
    const auto K = G.numberOfSpatialOrbitals();
    const auto orbital_space = OrbitalSpace::Implicit({{OccupationType::k_occupied, N_P}, {OccupationType::k_virtual, K - N_P}});  // N_P occupied (spatial) orbitals, K-N_P virtual (spatial) orbitals


    const auto& g = sq_hamiltonian.twoElectron().parameters();

    auto B_lambda = orbital_space.initializeRepresentableObjectFor<double>(OccupationType::k_occupied, OccupationType::k_virtual, OccupationType::k_occupied, OccupationType::k_virtual);  // initialize a mathematical object for the representaion of an occupied-virtual-occupied-virtual object (iabj)

    for (const auto& i : orbital_space.indices(OccupationType::k_occupied)) {
        for (const auto& a : orbital_space.indices(OccupationType::k_virtual)) {
            for (const auto& j : orbital_space.indices(OccupationType::k_occupied)) {
                for (const auto& b : orbital_space.indices(OccupationType::k_virtual)) {
                    double value {0.0};

                    value += lambda(i, b) * g(j, a, j, a) + lambda(j, a) * g(i, b, i, b);

                    if (i == j) {
                        value -= 2 * (lambda(i, b) * g(i, a, i, a) + lambda(i, a) * g(i, b, i, b));
                    }

                    if (a == b) {
                        value -= 2 * (lambda(j, a) * g(i, a, i, a) + lambda(i, a) * g(j, a, j, a));
                    }

                    if ((i == j) && (a == b)) {
                        value += 4 * lambda(i, a) * g(i, a, i, a);
                    }

                    B_lambda(i, a, j, b) = value;
                }
            }
        }
    }

    return B_lambda;
}


/**
 *  @param sq_hamiltonian                   the Hamiltonian expressed in an orthonormal orbital basis
 *  @param dipole_op                        the dipole integrals expressed in an orthonormal orbital basis
 *  @param x                                the first-order parameter response
 * 
 *  @return the Lagrangian multiplier response force (F_lambda)
 */
Matrix<double, Dynamic, 3> vAP1roGElectricalResponseSolver::calculateMultiplierResponseForce(const RSQHamiltonian<double>& sq_hamiltonian, const VectorRSQOneElectronOperator<double> dipole_op, const Matrix<double, Dynamic, 3>& x) const {

    const auto A_lambda = this->calculateExplicitMultiplierResponseForce(dipole_op);
    const auto B_lambda = this->calculateImplicitMultiplierResponseForceConstant(sq_hamiltonian);

    // Calculate F_lambda as F_lambda = A_lambda + B_lambda x
    return A_lambda + B_lambda.asMatrix() * x;
}


/**
 *  Solve the linear response equations for the Lagrangian multiplier response
 * 
 *  @param sq_hamiltonian               the Hamiltonian parameters
 *  @param dipole_op                    the dipole integrals expressed in an orthonormal orbital basis
 *  @param x                            the first-order parameter response
 * 
 *  @return the Lagrangian multiplier reponse y
 */
Matrix<double, Dynamic, 3> vAP1roGElectricalResponseSolver::calculateMultiplierResponse(const RSQHamiltonian<double>& sq_hamiltonian, const VectorRSQOneElectronOperator<double> dipole_op, const Matrix<double, Dynamic, 3>& x) const {

    const auto k_lambda = this->calculateMultiplierResponseConstant(sq_hamiltonian);
    const auto F_lambda = this->calculateMultiplierResponseForce(sq_hamiltonian, dipole_op, x);


    // This function is basically a wrapper around solving k_lambda y = -F_lambda
    auto environment = LinearEquationEnvironment<double>(k_lambda, -F_lambda);
    auto solver = LinearEquationSolver<double>::HouseholderQR();
    solver.perform(environment);

    const Matrix<double, Dynamic, 3> y = environment.x;
    return y;
}


}  // namespace GQCP
