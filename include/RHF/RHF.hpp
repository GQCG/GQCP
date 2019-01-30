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
#ifndef RHF_hpp
#define RHF_hpp


#include "HamiltonianParameters/HamiltonianParameters.hpp"
#include "RDM/OneRDM.hpp"

#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>


namespace GQCP {

/**
 *  A class that represents a converged solution to the RHF SCF equations
 */
class RHF {
private:
    double electronic_energy;
    Eigen::MatrixXd C;  // transformation matrix from the AO basis to the RHF MO basis
    Eigen::VectorXd orbital_energies;  // sorted in ascending energies


public:
    // CONSTRUCTORS
    /**
     *  Default constructor setting everything to zero
     */
    RHF();  // need default constructor

    /**
     *  Constructor based on given converged solutions of the RHF SCF equations
     *
     *  @param electronic_energy    the converged RHF electronic energy
     *  @param C                    the coefficient matrix, i.e. the transformation matrix from the AO basis to the RHF MO basis
     *  @param orbital_energies     the RHF MO energies
     */
    RHF(double electronic_energy, const Eigen::MatrixXd& C, const Eigen::VectorXd& orbital_energies);


    // GETTERS
    double get_electronic_energy() const { return this->electronic_energy; }
    const Eigen::MatrixXd& get_C() const { return this->C; }
    const Eigen::VectorXd& get_orbital_energies() const { return this->orbital_energies; }
    double get_orbital_energies(size_t index) const { return this->orbital_energies(index); }
};


/*
 *  HELPER METHODS
 */

/**
 *  @param K    the number of spatial orbitals
 *  @param N    the number of electrons
 *
 *  @return the RHF 1-RDM expressed in an orthonormal basis
 */
OneRDM calculateRHF1RDM(size_t K, size_t N);

/**
 *  @param C    the coefficient matrix, specifying the transformation to the AO basis
 *  @param N    the number of electrons
 *
 *  @return the RHF 1-RDM expressed in the AO basis
 */
Eigen::MatrixXd calculateRHFAO1RDM(const Eigen::MatrixXd& C, size_t N);

/**
 *  Calculate the RHF Fock matrix F = H_core + G, in which G is a contraction of the density matrix and the two-electron integrals
 *
 *  @param D_AO     the RHF density matrix in AO basis
 *  @param ham_par  The Hamiltonian parameters in AO basis
 *
 *  @return the RHF Fock matrix expressed in the AO basis
 */
Eigen::MatrixXd calculateRHFAOFockMatrix(const Eigen::MatrixXd& D_AO, HamiltonianParameters ham_par);

/**
 *  @param D_AO         the RHF density matrix in AO basis
 *  @param H_core_AO    the core Hamiltonian parameters in AO basis
 *  @param F_AO         the Fock matrix in AO basis
 *
 *  @return the RHF electronic energy
 */
double calculateRHFElectronicEnergy(const Eigen::MatrixXd& D_AO, const Eigen::MatrixXd& H_core_AO, const Eigen::MatrixXd& F_AO);

/**
 *  @param N    the number of electrons
 *
 *  @return the RHF HOMO index
 */
size_t RHFHOMOIndex(size_t N);

/**
 *  @param K    the number of spatial orbitals
 *  @param N    the number of electrons
 *
 *  @return the RHF LUMO index
 */
size_t RHFLUMOIndex(size_t K, size_t N);


}  // namespace GQCP


#endif /* RHF_hpp */
