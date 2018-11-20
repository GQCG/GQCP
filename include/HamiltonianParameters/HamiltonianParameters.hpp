// This file is part of GQCG-gqcp.
// 
// Copyright (C) 2017-2018  the GQCG developers
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
#ifndef GQCP_HAMILTONIANPARAMETERS_HPP
#define GQCP_HAMILTONIANPARAMETERS_HPP

#include "common.hpp"
#include "HamiltonianParameters/BaseHamiltonianParameters.hpp"
#include "Operator/OneElectronOperator.hpp"
#include "Operator/TwoElectronOperator.hpp"
#include "JacobiRotationParameters.hpp"
#include "RDM/TwoRDM.hpp"
#include "RDM/OneRDM.hpp"

#include <Eigen/Dense>



namespace GQCP {


/**
 *  A class for representing Hamiltonian parameters, i.e. the one- and two-electron integrals in the second-quantized expression of the Hamiltonian
 *
 *  This class can be used for restricted calculations, i.e. the alpha and beta integrals are equal
 */
class HamiltonianParameters : public BaseHamiltonianParameters {
private:
    size_t K;  // the number of spatial orbitals

    OneElectronOperator S;  // overlap

    OneElectronOperator h;  // one-electron interactions (i.e. the core Hamiltonian)
    TwoElectronOperator g;  // two-electron interactions

    Eigen::MatrixXd C;  // total transformation matrix between the current (restricted) molecular orbitals and the atomic orbitals


public:
    // CONSTRUCTORS
    /**
     *  @param ao_basis     the initial AO basis
     *  @param S            the overlap integrals
     *  @param h            the one-electron integrals H_core
     *  @param g            the two-electron integrals
     *  @param C            a transformation matrix between the current molecular orbitals and the atomic orbitals
     *  @param scalar       the scalar interaction term
     */
    HamiltonianParameters(std::shared_ptr<GQCP::AOBasis> ao_basis, const GQCP::OneElectronOperator& S, const GQCP::OneElectronOperator& h, const GQCP::TwoElectronOperator& g, const Eigen::MatrixXd& C, double scalar=0.0);


    /**
     *  A constructor that transforms the current Hamiltonian parameters with a transformation matrix
     *
     *  @param ham_par      the current Hamiltonian parameters
     *  @param C            the transformation matrix to be applied to the current Hamiltonian parameters
     */
    HamiltonianParameters(const GQCP::HamiltonianParameters& ham_par, const Eigen::MatrixXd& C);


    // DESTRUCTORS
    ~HamiltonianParameters() override =default;

    
    // GETTERS
    const GQCP::OneElectronOperator& get_S() const { return this->S; }
    const GQCP::OneElectronOperator& get_h() const { return this->h; }
    const GQCP::TwoElectronOperator& get_g() const { return this->g; }
    const Eigen::MatrixXd& get_C() const { return this->C; }
    size_t get_K() const { return this->K; }

    
    // PUBLIC METHODS
    /**
     *  In-place transform the matrix representations of Hamiltonian parameters
     *
     *  @param T    the transformation matrix between the old and the new orbital basis, it is used as
     *      b' = b T ,
     *   in which the basis functions are collected as elements of a row vector b
     *
     *  Furthermore
     *      - the overlap matrix S now gives the overlap matrix in the new molecular orbital basis
     *      - the coefficient matrix C is updated to reflect the total transformation between the new molecular orbital basis and the initial atomic orbitals
     */
    void transform(const Eigen::MatrixXd& T);

    /**
     *  In-place rotate the matrix representations of the Hamiltonian parameters
     *
     *  @param U     the unitary transformation (i.e. rotation) matrix, see transform() for how the transformation matrix between the two bases should be represented
     *
     *  Furthermore, the coefficient matrix C is updated to reflect the total transformation between the new molecular orbital basis and the initial atomic orbitals
     */
    void rotate(const Eigen::MatrixXd& U);

    /**
     *  Using a random rotation matrix, transform the matrix representations of the Hamiltonian parameters
     */
    void randomRotate();

    /**
     *  In-place rotate the matrix representations of the Hamiltonian parameters using a unitary Jacobi rotation matrix constructed from the Jacobi rotation parameters
     *
     *  @param jacobi_rotation_parameters       the Jacobi rotation parameters (p, q, angle) that are used to specify a Jacobi rotation: we use the (cos, sin, -sin, cos) definition for the Jacobi rotation matrix. See transform() for how the transformation matrix between the two bases should be represented
     *
     *  Furthermore the coefficient matrix C is updated to reflect the total transformation between the new molecular orbital basis and the initial atomic orbitals
     */
    void rotate(const GQCP::JacobiRotationParameters& jacobi_rotation_parameters);

    /**
     *  Transform the HamiltonianParameters to the Löwdin basis (i.e. T = S^{-1/2})
     */
    void LowdinOrthonormalize();

    /**
     *  @param D      the 1-RDM
     *  @param d      the 2-RDM
     *
     *  @return the generalized Fock matrix
     */
    GQCP::OneElectronOperator calculateGeneralizedFockMatrix(const GQCP::OneRDM& D, const GQCP::TwoRDM& d) const;

    /**
     *  @param D      the 1-RDM
     *  @param d      the 2-RDM
     *
     *  @return the super-generalized Fock matrix
     */
    GQCP::TwoElectronOperator calculateSuperGeneralizedFockMatrix(const GQCP::OneRDM& D, const GQCP::TwoRDM& d) const;

    /**
     *  @param N_P      the number of electron pairs
     *
     *  @return the Edmiston-Ruedenberg localization index g(i,i,i,i)
     */
    double calculateEdmistonRuedenbergLocalizationIndex(size_t N_P) const;

    /**
     *  Constrain the Hamiltonian parameters according to the convention: - lambda * constraint
     *
     *  @param one_op   the one-electron operator used as a constraint
     *  @param two_op   the two-electron operator used as a constraint
     *  @param lambda   Lagrangian multiplier for the constraint
     *
     *  @return a copy of the constrained Hamiltonian parameters
     */
    HamiltonianParameters constrain(const GQCP::OneElectronOperator& one_op, const GQCP::TwoElectronOperator& two_op, double lambda) const;

    /**
     *  Constrain the Hamiltonian parameters according to the convention: - lambda * constraint
     *
     *  @param one_op   the one-electron operator used as a constraint
     *  @param lambda   Lagrangian multiplier for the constraint
     *
     *  @return a copy of the constrained Hamiltonian parameters
     */
    HamiltonianParameters constrain(const GQCP::OneElectronOperator& one_op, double lambda) const;

    /**
     *  Constrain the Hamiltonian parameters according to the convention: - lambda * constraint
     *
     *  @param two_op   the two-electron operator used as a constraint
     *  @param lambda   Lagrangian multiplier for the constraint
     *
     *  @return a copy of the constrained Hamiltonian parameters
     */
    HamiltonianParameters constrain(const GQCP::TwoElectronOperator& two_op, double lambda) const;

    /**
     *  @param ao_list     indices of the AOs used for the Mulliken populations
     *
     *  @return the Mulliken operator for a set of AOs
     */
    OneElectronOperator calculateMullikenOperator(const Vectoru& ao_list);
};


}  // namespace GQCP


#endif  // GQCP_HAMILTONIANPARAMETERS_HPP
