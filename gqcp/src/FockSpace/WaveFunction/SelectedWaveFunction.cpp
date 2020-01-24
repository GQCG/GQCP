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
#include "FockSpace/WaveFunction/SelectedWaveFunction.hpp"


namespace GQCP {


/*
 * CONSTRUCTORS
 */

/**
 *  Construct a selected wave function
 *
 *  @param fock_space           the selected Fock space in which the wave function 'lives'
 *  @param coefficients         the expansion coefficients
 */
SelectedWaveFunction::SelectedWaveFunction(const SelectedFockSpace& fock_space, const VectorX<double>& coefficients)
    fock_space (fock_space),
    coefficients (coefficients)
{
    if (std::abs(this->coefficients.norm() - 1.0) > 1.0e-12) {  // normalize the coefficients if they aren't
        this->coefficients.normalize();
    }
}



/**
 *  Construct a selected wave function by selecting the configurations with the highest contribution from a different wave function
 *
 *  @param wave_function                           the wave function to be selected from
 *  @param number_of_highest_contributions         the amount of selections made
 */
SelectedWaveFunction::SelectedWaveFunction(const WaveFunction& wave_function, size_t number_of_highest_contributions) {

    struct AddressCoefficientPair 
    {
        double coeff;
        size_t index;
    };

    struct AddressCoefficientPairComparer 
    { 
        int operator() (const Pair& p1, const Pair& p2) 
        { 
            return std::abs(p1.coeff) >std::abs(p2.coeff); 
        } 
    }; 

    // Initialize a "min heap" datastructure for which the top value is the lowest value
   std::priority_queue <AddressCoefficientPair, std::vector<AddressCoefficientPair>, AddressCoefficientPairComparer > min_heap; 
   
   std::ofstream outfile ("berylium2_2x.txt");
   outfile<<std::setprecision(16);
   GQCP::Nucleus Be (GQCP::elements::elementToAtomicNumber("Be"), 0, 0, 0);
   GQCP::Molecule Be_mol ({Be}, -2);
   GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis (Be_mol, "aug-cc-pvdz");

   auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, Be_mol);  // in an AO basis
   auto K = sq_hamiltonian.dimension();
   outfile << "K:"<< K  << std::endl;
   // Create a plain RHF SCF solver and solve the SCF equations
   GQCP::DIISRHFSCFSolver diis_scf_solver (sq_hamiltonian, spinor_basis, Be_mol);
   diis_scf_solver.solve();
   auto rhf = diis_scf_solver.get_solution();
   outfile << "RHF Orbital energies:"<< std::endl << rhf.get_orbital_energies() << std::endl;

   // Transform the Hamiltonian to the RHF orbital basis
   GQCP::basisTransform(spinor_basis, sq_hamiltonian, rhf.get_C());

   GQCP::ProductFockSpace fock_space (K, Be_mol.numberOfElectrons()/2, Be_mol.numberOfElectrons()/2);  
   const auto& fock_space_alpha = fock_space.get_fock_space_alpha();
   const auto& fock_space_beta = fock_space.get_fock_space_beta();

   // Create the FCI module
   GQCP::FCI fci (fock_space);
   GQCP::CISolver ci_solver (fci, sq_hamiltonian);

   // Solve Davidson
   GQCP::VectorX<double> initial_g = fock_space.HartreeFockExpansion();
   GQCP::DavidsonSolverOptions davidson_solver_options (initial_g);
   ci_solver.solve(davidson_solver_options);

   // Retrieve the eigenvalues
   auto fci_davidson_eigenvalue = ci_solver.get_eigenpair().get_eigenvalue();
   auto energy = fci_davidson_eigenvalue + GQCP::Operator::NuclearRepulsion(Be_mol).value();
   auto wave = ci_solver.makeWavefunction();
   GQCP::RDMCalculator rdm_calculator(wave);
   GQCP::TwoRDMs<double> t = rdm_calculator.calculate2RDMs();
   GQCP::OneRDMs<double> o = rdm_calculator.calculate1RDMs();

   auto s2 = GQCP::calculateSpinSquared(o, t);

   outfile << "Eigenvalue:"<< energy << std::endl;
   outfile << "<S^2>:"<< s2 << std::endl;

   Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> saes(o.one_rdm);
   auto U = saes.eigenvectors(); 
   outfile << "D_a:"<< std::endl << o.one_rdm_aa << std::endl;
   Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> saes2(o.one_rdm_aa);
   outfile << "Y_a:"<< std::endl << saes2.eigenvalues() << std::endl;
   outfile << "D_b:"<< std::endl << o.one_rdm_bb << std::endl;
   Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> saes3(o.one_rdm_bb);
   outfile << "Y_b:"<< std::endl << saes3.eigenvalues() << std::endl;

   outfile << "D_d:"<< std::endl << o.spinDensityRDM() << std::endl;
   Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> saes4( o.spinDensityRDM() );
   outfile << "Y_d:"<< std::endl << saes4.eigenvalues() << std::endl;

   outfile << "D:"<< std::endl << o.one_rdm << std::endl;
   outfile << "Y:"<< std::endl << saes.eigenvalues() << std::endl;

   outfile << "RHF coefficients:"<< std::endl << rhf.get_C()<< std::endl;
   outfile << "T coefficients:"<< std::endl << rhf.get_C() * U << std::endl;

   GQCP::basisTransform(spinor_basis, sq_hamiltonian, GQCP::TransformationMatrix<double>(U));

   const auto& total_one_op = sq_hamiltonian.core();
   const auto& total_two_op = sq_hamiltonian.twoElectron();
   size_t K2 = 9;
   size_t K3 = 18;
   const auto& chipped_one = total_one_op.parameters().block(K-K2, K-K2, K2, K2);
   const auto& chipped_one2 = total_one_op.parameters().block(K-K3, K-K3, K3, K3);
   const auto& chipped_two = GQCP::Tensor<double, 4>::FromBlock(total_two_op.parameters(),  K-K2, K-K2, K-K2, K-K2);
   const auto& chipped_two2 =  GQCP::Tensor<double, 4>::FromBlock(total_two_op.parameters(),  K-K3, K-K3, K-K3, K-K3);

   // SUB 
   GQCP::SQHamiltonian<double> nat_ham ({GQCP::ScalarSQOneElectronOperator<double>({chipped_one})}, {GQCP::ScalarSQTwoElectronOperator<double>({chipped_two})});
   GQCP::SQHamiltonian<double> nat_ham2 ({GQCP::ScalarSQOneElectronOperator<double>({chipped_one2})}, {GQCP::ScalarSQTwoElectronOperator<double>({chipped_two2})});
   GQCP::ProductFockSpace fock_space2 (K2, Be_mol.numberOfElectrons()/2, Be_mol.numberOfElectrons()/2);  
   GQCP::ProductFockSpace fock_space3 (K3, Be_mol.numberOfElectrons()/2, Be_mol.numberOfElectrons()/2);  
   const auto& fock_space_alpha2 = fock_space2.get_fock_space_alpha();
   const auto& fock_space_alpha3 = fock_space3.get_fock_space_alpha();
   const auto& fock_space_beta2 = fock_space2.get_fock_space_beta();
   const auto& fock_space_beta3 = fock_space3.get_fock_space_beta();

   // Create the FCI module
   GQCP::FCI fci2 (fock_space2);
   GQCP::FCI fci3 (fock_space3);
   GQCP::CISolver ci_solver2 (fci2, nat_ham);

   // Solve Davidson
   GQCP::DenseSolverOptions dense;
   ci_solver2.solve(dense);
   auto wave2 = ci_solver2.makeWavefunction();
   const auto& coeff = wave2.get_coefficients();

   GQCP::VectorX<double> guess_i = GQCP::VectorX<double>::Zero(fock_space.get_dimension());

   for (size_t i = 0; i < fock_space2.get_dimension(); i++) {
      size_t alpha_index = i/fock_space_alpha2.get_dimension();
      size_t beta_index = i - alpha_index * fock_space_alpha2.get_dimension();
      const GQCP::ONV alpha = fock_space_alpha2.makeONV(alpha_index);
      const GQCP::ONV beta = fock_space_beta2.makeONV(beta_index);

      size_t address = fock_space_alpha.getAddress(alpha.get_unsigned_representation()  << (K - K2) ) * fock_space_alpha.get_dimension() + fock_space_beta.getAddress(beta.get_unsigned_representation() << (K - K2));
      guess_i(address) = coeff(i);

   }

   /*
   GQCP::CISolver ci_solver3 (fci3, nat_ham2);
   davidson_solver_options.X_0 = guess_i;
   ci_solver3.solve(davidson_solver_options);
   auto wave3 = ci_solver3.makeWavefunction();
   const auto coefficients3 = wave3.get_coefficients();

   GQCP::VectorX<double> guess = GQCP::VectorX<double>::Zero(fock_space.get_dimension());

   for (size_t i = 0; i < fock_space2.get_dimension(); i++) {
      size_t alpha_index = i/fock_space_alpha3.get_dimension();
      size_t beta_index = i - alpha_index * fock_space_alpha3.get_dimension();
      const GQCP::ONV alpha = fock_space_alpha3.makeONV(alpha_index);
      const GQCP::ONV beta = fock_space_beta3.makeONV(beta_index);

      size_t address = fock_space_alpha.getAddress(alpha.get_unsigned_representation() << (K-K3)) * fock_space_alpha.get_dimension() + fock_space_beta.getAddress(beta.get_unsigned_representation()  << (K-K3));
      guess(address) = coefficients3(i);

   }
   */

  davidson_solver_options.maximum_subspace_dimension = 5;

   GQCP::CISolver ci_solver4 (fci, sq_hamiltonian);
   davidson_solver_options.X_0 = guess_i;
   ci_solver4.solve(davidson_solver_options);
   auto wave4 = ci_solver4.makeWavefunction();
   const auto coefficients3 = wave4.get_coefficients();
   for (size_t i = 0; i < fock_space.get_dimension(); i++) {
      pq.push(Pair {coefficients3(i), i});
      if (pq.size() > 1000) {
         pq.pop();
      }
   }

   while (!pq.empty()) {
      const Pair& x = pq.top();
      size_t alpha_index = x.index/fock_space_alpha.get_dimension();
      size_t beta_index = x.index - alpha_index * fock_space_alpha.get_dimension();
      const GQCP::ONV alpha = fock_space_alpha.makeONV(alpha_index);
      const GQCP::ONV beta = fock_space_beta.makeONV(beta_index);
      outfile << alpha.asString() << " | " << beta.asString() << " : "<< x.coeff << " sen:  " << alpha.countNumberOfDifferences(beta) << std::endl;
      pq.pop();
   }

   outfile.close();
   return 0;


}


// PUBLIC METHODS
/**
 *  Print this wave function to an output in the GAMESS-US Expansion format
 *
 *  @param output_stream        the stream used for outputting
 */
void SelectedWaveFunction::printGamessExpansion(std::ostream& output_stream = std::cout) const;


}  // namespace GQCP
