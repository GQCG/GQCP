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
#pragma once


#include "Basis/TransformationMatrix.hpp"
#include "Mathematical/Representation/BlockRankFourTensor.hpp"
#include "Mathematical/Representation/Tensor.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "Processing/RDM/OneRDM.hpp"


namespace GQCP {


/**
 *  A class that represents a converged solution to the RHF SCF equations
 */
class RHF {
private:
    double electronic_energy;
    TransformationMatrix<double> C;  // transformation matrix from the AO basis to the RHF MO basis
    VectorX<double> orbital_energies;  // sorted in ascending energies


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
    RHF(double electronic_energy, const TransformationMatrix<double>& C, const VectorX<double>& orbital_energies);


    // GETTERS
    double get_electronic_energy() const { return this->electronic_energy; }
    const TransformationMatrix<double>& get_C() const { return this->C; }
    const VectorX<double>& get_orbital_energies() const { return this->orbital_energies; }
    double get_orbital_energies(size_t index) const { return this->orbital_energies(index); }
};


/*
 *  HELPER METHODS
 */

/**
 *  Calculate the RHF Fock matrix F = H_core + G, in which G is a contraction of the density matrix and the two-electron integrals
 *
 *  @param D_AO                 the RHF density matrix in AO basis
 *  @param sq_hamiltonian       the Hamiltonian expressed in an AO basis
 *
 *  @return the RHF Fock matrix expressed in the AO basis
 */
ScalarSQOneElectronOperator<double> calculateRHFAOFockMatrix(const OneRDM<double>& D_AO, const SQHamiltonian<double>& sq_hamiltonian);

/**
 *  @param D_AO         the RHF density matrix in AO basis
 *  @param H_core_AO    the core Hamiltonian expressed in an AO basis
 *  @param F_AO         the Fock matrix in AO basis
 *
 *  @return the RHF electronic energy
 */
double calculateRHFElectronicEnergy(const OneRDM<double>& D_AO, const ScalarSQOneElectronOperator<double>& H_core_AO, const ScalarSQOneElectronOperator<double>& F_AO);

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

/**
 *  Specialize the orbital Hessian for RHF
 * 
 *  @param sq_hamiltonian       the Hamiltonian expressed in an orthonormal basis
 *  @param N_P                  the number of electron pairs
 * 
 *  @return the RHF orbital Hessian as a BlockRankFourTensor, i.e. an object with a suitable operator() implemented
 */
BlockRankFourTensor<double> calculateRHFOrbitalHessianTensor(const SQHamiltonian<double>& sq_hamiltonian, const size_t N_P);

/**
 *  @param sq_hamiltonian       the Hamiltonian expressed in an orthonormal basis
 *  @param N_P                  the number of electron pairs
 *  @param a                    the first virtual orbital index
 *  @param i                    the first occupied orbital index
 *  @param b                    the second virtual orbital index
 *  @param j                    the second occupied orbital index
 * 
 *  @return an element of the RHF orbital Hessian
 */
double calculateRHFOrbitalHessianElement(const SQHamiltonian<double>& sq_hamiltonian, const size_t N_P, const size_t a, const size_t i, const size_t b, const size_t j);


}  // namespace GQCP
