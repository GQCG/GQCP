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


#include "Basis/SpinorBasis/RSpinorBasis.hpp"
#include "CISolver/CISolver.hpp"
#include "FockSpace/ProductFockSpace.hpp"
#include "HamiltonianBuilder/FCI.hpp"
#include "HamiltonianBuilder/FrozenCoreFCI.hpp"
#include "HamiltonianParameters/AtomicDecompositionParameters.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "Operator/SecondQuantized/USQHamiltonian.hpp"
#include "RDM/RDMCalculator.hpp"


namespace GQCP {
namespace QCMethod {


/**
 *  A class that solves the FCI Hamiltonian given a perturbation in the form of a langragian multiplier and the Mulliken operator for a pre-specified set of basis functions
 *  Additionally an atomic Sz perturbation can be applied as wel. 
 */
class MullikenConstrainedFCI {
private:
    double solve_time;
    std::vector<size_t> basis_targets;
    Molecule molecule;
    RSpinorBasis<double, GTOShell> spinor_basis;
    USpinorBasis<double, GTOShell> uspinor_basis;
    SQHamiltonian<double> sq_hamiltonian;
    USQHamiltonian<double> usq_hamiltonian;
    std::string basis_set;  // the basisset that should be used
    FrozenProductFockSpace fock_space = FrozenProductFockSpace(0, 0, 0, 0); // Default
    FrozenCoreFCI fci = FrozenCoreFCI(FrozenProductFockSpace(0, 0, 0, 0)); 
    ScalarSQOneElectronOperator<double> mulliken_operator;
    ScalarSQOneElectronOperator<double> sq_sz_operator;
    AtomicDecompositionParameters adp = AtomicDecompositionParameters();
    RDMCalculator rdm_calculator = RDMCalculator();

    bool are_solutions_available = false;

    // Molecular solutions
    std::vector<double> energy;
    std::vector<double> population; 
    std::vector<double> lambda;
    std::vector<double> lambda_sz;
    std::vector<double> entropy;
    std::vector<double> sz;


    // Decomposed solutions
    std::vector<double> A_fragment_energy;
    std::vector<double> A_fragment_self_energy;
    std::vector<double> B_fragment_energy;
    std::vector<double> B_fragment_self_energy;
    std::vector<double> interaction_energy;

    // Eigenvectors
    std::vector<VectorX<double>> eigenvector;

    // PRIVATE METHODS
    /**
     *  Store the solutions from a solve
     *  
     *  @param eigenpairs           the eigenpairs from the CI solver
     *  @param multiplier           the Lagrangian multiplier associated with the solution
     *  @param sz_multiplier        a given multiplier for the atomic Sz constraint
     */ 
    void parseSolution(const std::vector<Eigenpair>& eigenpairs, const double multiplier, const double sz_multiplier = 0);

    /**
     *  Throws an error if no solution is available
     *  
     *  @param function_name            name of the function that should throw the error
     */
    void checkAvailableSolutions(const std::string& function_name) const;

    /**
     *  Throws an error if the molecule is not diatomic
     *  
     *  @param function_name            name of the function that should throw the error
     */
    void checkDiatomicMolecule(const std::string& function_name) const;

public:

    // CONSTRUCTORS
    /**
     *  @param molecule                 the molecule that will be solved for
     *  @param basis_set                the basisset that should be used
     *  @param basis_targets            the targeted basis functions for the constraint
     *  @param frozencores              the amount of frozen cores for the FCI calculation
     */
    MullikenConstrainedFCI(const Molecule& molecule, const std::string& basis_set, const std::vector<size_t>& basis_targets, const size_t frozencores = 0);


    // PUBLIC METHODS
    /**
     *  Solve the eigenvalue problem for a multiplier with the davidson algorithm
     *  
     *  @param multiplier           a given multiplier for the Mulliken constraint
     *  @param guess                supply a davidson guess
     *  @param sz_multiplier        a given multiplier for the atomic Sz constraint
     */
    void solveMullikenDavidson(const double multiplier, const VectorX<double>& guess, const double sz_multiplier = 0);

    /**
     *  Solve the eigenvalue problem for a multiplier with the davidson algorithm, davidson guess will be the previously stored solution
     *  if none is available the Hartree Fock expansion will be used instead
     *  
     *  @param multiplier           a given multiplier
     *  @param sz_multiplier        a given multiplier for the atomic Sz constraint
     */
    void solveMullikenDavidson(const double multiplier, const double sz_multiplier = 0);

    /**
     *  Solve the eigenvalue problem for a the next multiplier dense
     * 
     *  @param multiplier           a given multiplier
     *  @param nos                  the number of eigenpairs or "states" that should be stored for each multiplier``
     *  @param sz_multiplier        a given multiplier for the atomic Sz constraint
     */
    void solveMullikenDense(const double multiplier, const size_t nos, const double sz_multiplier = 0);

    /**
     *  @param index                refers to the index of the number of requested states 
     * 
     *  @return a property of the last solve
     */
    double get_energy(const size_t index = 0) const { this->checkAvailableSolutions("get_energy"); return this->energy[index]; };
    double get_population(const size_t index = 0) const { this->checkAvailableSolutions("get_population"); return this->population[index]; };
    double get_lambda(const size_t index = 0) const { this->checkAvailableSolutions("get_lambda"); return this->lambda[index]; };
    double get_lambda_sz(const size_t index = 0) const { this->checkAvailableSolutions("get_lambda_sz"); return this->lambda_sz[index]; };
    double get_entropy(const size_t index = 0) const { this->checkAvailableSolutions("get_entropy"); return this->entropy[index]; };
    double get_sz(const size_t index = 0) const { this->checkAvailableSolutions("get_sz"); return this->sz[index]; };
    double get_A_fragment_energy(const size_t index = 0) const { this->checkAvailableSolutions("get_A_fragment_energy"); this->checkDiatomicMolecule("get_A_fragment_energy"); return this->A_fragment_energy[index]; };
    double get_A_fragment_self_energy(const size_t index = 0) const { this->checkAvailableSolutions("get_A_fragment_self_energy"); this->checkDiatomicMolecule("get_A_fragment_self_energy"); return this->A_fragment_self_energy[index]; };
    double get_B_fragment_energy(const size_t index = 0) const { this->checkAvailableSolutions("get_B_fragment_energy"); this->checkDiatomicMolecule("get_B_fragment_energy"); return this->B_fragment_energy[index]; };
    double get_B_fragment_self_energy(const size_t index = 0) const { this->checkAvailableSolutions("get_B_fragment_self_energy"); this->checkDiatomicMolecule("get_B_fragment_self_energy"); return this->B_fragment_self_energy[index]; };
    double get_interaction_energy(const size_t index = 0) const { this->checkAvailableSolutions("get_interaction_energy"); this->checkDiatomicMolecule("get_interaction_energy"); return this->interaction_energy[index]; };

    const VectorX<double>& get_eigenvector(const size_t index = 0) const { this->checkAvailableSolutions("get_eigenvector"); return this->eigenvector[index]; };

    double get_solve_time() const { this->checkAvailableSolutions("get_solve_time"); return solve_time; };
    
    /**
     *  @param index             refers to the index of the number of requested states 
     * 
     *  @return all properties in vector that contains:
     *      energy, population (on the selected basis functions), Sz, lambda (or the multiplier for the Mulliken constraint), lambda_sz (for atomic Sz), entropy
     *      if diatomic we additionally find: A_fragment_energy, A_fragment_self_energy, B_fragment_energy, B_fragment_self_energy and interaction_energy in that order.
     */
    std::vector<double> all_properties(const size_t index = 0) const; 
};


}  // namespace QCMethod
}  // namespace GQCP
