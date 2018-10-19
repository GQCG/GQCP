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


#include "HamiltonianParameters/BaseHamiltonianParameters.hpp"
#include "Operator/OneElectronOperator.hpp"
#include "Operator/TwoElectronOperator.hpp"
#include "JacobiRotationParameters.hpp"
#include "RDM/TwoRDM.hpp"
#include "RDM/OneRDM.hpp"

#include <Eigen/Dense>



namespace GQCP {


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
     *  Constructor based on a given @param ao_basis, overlap @param S, one-electron operator @param h, two-electron
     *  operator @param g and a transformation matrix between the current molecular orbitals and the atomic orbitals
     *  @param C
     */
    HamiltonianParameters(std::shared_ptr<GQCP::AOBasis> ao_basis, const GQCP::OneElectronOperator& S, const GQCP::OneElectronOperator& h, const GQCP::TwoElectronOperator& g, const Eigen::MatrixXd& C);


    /**
     *  Constructor based on given Hamiltonian parameters @param ham_par and a transformation matrix @param C.
     *
     *  If the initial Hamiltonian parameters @param ham_par are expressed in the basis B, the constructed instance represents the Hamiltonian parameters in the transformed basis B'. The basis transformation between B and B' is given by the transformation matrix @param C.
     */
    HamiltonianParameters(const GQCP::HamiltonianParameters& ham_par, const Eigen::MatrixXd& C);


    // DESTRUCTORS
    ~HamiltonianParameters() override =default;

    
    // GETTERS
    GQCP::OneElectronOperator get_h() const { return this->h; }
    GQCP::TwoElectronOperator get_g() const { return this->g; }
    size_t get_K() const { return this->K; }

    
    // PUBLIC METHODS
    /**
     *  Given a transformation matrix @param T that links the new molecular orbital basis to the old molecular orbital basis,
     *  in the sense that
     *       b' = b T ,
     *  in which the molecular orbitals are collected as elements of a row vector b, transform
     *      - the one-electron interaction operator (i.e. the core Hamiltonian)
     *      - the two-electron interaction operator
     *
     *  Furthermore
     *      - @member S now gives the overlap matrix in the new molecular orbital basis
     *      - @member C is updated to reflect the total transformation between the new molecular orbital basis and the initial atomic orbitals
     */
    void transform(const Eigen::MatrixXd& T);

    /**
     *  Given a unitary rotation matrix @param U that links the new molecular orbital basis to the old molecular orbital basis,
     *  in the sense that
     *       b' = b U ,
     *  in which the molecular orbitals are collected as elements of a row vector b, transform
     *      - the one-electron interaction operator (i.e. the core Hamiltonian)
     *      - the two-electron interaction operator
     *
     *  Furthermore, @member C is updated to reflect the total transformation between the new molecular orbital basis and the initial atomic orbitals
     */
    void rotate(const Eigen::MatrixXd& U);

    /**
     *  Given @param jacobi_rotation_parameters that represent a unitary rotation matrix @param U (using a (cos, sin, -sin, cos) definition for the Jacobi rotation matrix) that links the new molecular orbital basis to the old molecular orbital basis,
     *  in the sense that
     *       b' = b U ,
     *  in which the molecular orbitals are collected as elements of a row vector b, transform
     *      - the one-electron interaction operator (i.e. the core Hamiltonian)
     *      - the two-electron interaction operator
     *
     *  Furthermore @member C is updated to reflect the total transformation between the new molecular orbital basis and the initial atomic orbitals
     */
    void rotate(const GQCP::JacobiRotationParameters& jacobi_rotation_parameters);

    /**
     *  Given @param one_rdm and @param two_rdm
     *  @return the energy as a result of the contraction of the 1- and 2-RDMs with the one- and two-electron integrals
     */
    double calculateEnergy(const GQCP::OneRDM& one_rdm, const GQCP::TwoRDM& two_rdm) const;

    /**
     *  Given a @param one_rdm and a @param two_rdm, @return the generalized Fock matrix F as a OneElectronOperator
     */
    GQCP::OneElectronOperator calculateGeneralizedFockMatrix(const GQCP::OneRDM& one_rdm, const GQCP::TwoRDM& two_rdm) const;

    /**
     *  Given a @param D: 1-RDM and a @param d: 2-RDM, @return the super-generalized Fock matrix W as a TwoElectronOperator
     */
    GQCP::TwoElectronOperator calculateSuperGeneralizedFockMatrix(const GQCP::OneRDM& one_rdm, const GQCP::TwoRDM& two_rdm) const;


    // FRIEND FUNCTIONS
    friend Eigen::MatrixXd calculateRHFAOFockMatrix(const Eigen::MatrixXd& D_AO, GQCP::HamiltonianParameters ham_par);
    friend double calculateRMP2EnergyCorrection(const GQCP::HamiltonianParameters& ham_par);

    // FRIEND CLASSES
    friend class RHFSCFSolver;
    friend class DIISRHFSCFSolver;
    friend class AP1roGPSESolver;
    friend class AP1roGJacobiOrbitalOptimizer;
};


}  // namespace GQCP


#endif  // GQCP_HAMILTONIANPARAMETERS_HPP
