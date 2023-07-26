# Import Statements.
import gqcpy
import numpy as np
import numpy.random as rand
import copy
import GQCC

from prototypes.helper import *

from GQCC.Operators.Operator import UnresolvedSpinorBasisOperator
from GQCC.Operators.Operator import ResolvedSpinorBasisOperator

from GQCC.Optimization.UnresolvedOptimizer import UnresolvedOptimizer


# Operator classes.
class RealSSquaredGOperator(UnresolvedSpinorBasisOperator):
    
    def __init__(self, basis):
        # Quantize the operator from the basis and store it.
        Ka = basis.numberOfCoefficients(gqcpy.Spin.alpha)
        
        S = basis.quantize(gqcpy.OverlapOperator()).parameters()
        S_aa = S[:Ka, :Ka]

        Sz = basis.quantize(gqcpy.ElectronicSpin_zOperator()).parameters()
        self.Sz = gqcpy.ScalarGSQOneElectronOperator_d(Sz)

        Sz2 = 0.25 * np.vstack((np.hstack((S_aa, np.zeros_like(S_aa))), np.hstack((np.zeros_like(S_aa), S_aa))))

        S_plus = np.vstack((np.hstack((np.zeros_like(S_aa), S_aa)), np.hstack((np.zeros_like(S_aa), np.zeros_like(S_aa)))))
        S_min = np.vstack((np.hstack((np.zeros_like(S_aa), np.zeros_like(S_aa))), np.hstack((S_aa, np.zeros_like(S_aa)))))
        S_min_plus = np.vstack((np.hstack((np.zeros_like(S_aa), np.zeros_like(S_aa))), np.hstack((np.zeros_like(S_aa), S_aa))))

        spin = Sz2 + Sz + S_min_plus
        tensor = np.einsum('pq,rs->pqrs', Sz, Sz) + np.einsum('pq,rs->pqrs', S_min, S_plus)

        self.matrix_operator = gqcpy.ScalarGSQOneElectronOperator_d(spin)
        self.tensor_operator = gqcpy.ScalarGSQTwoElectronOperator_d(tensor*2)

        
        # For later uses, the operator object contains the "name" of the resulting expectation values.
        self.operator_result = "S^2 value"


class RealSSquaredUOperator(ResolvedSpinorBasisOperator):

    def __init__(self, basis):
        # Quantize the operator from the basis and store it.
        # Overlaps.
        S_aa = basis.quantize(gqcpy.OverlapOperator()).alpha.parameters()
        S_bb = basis.quantize(gqcpy.OverlapOperator()).beta.parameters()

        # S_z operator.
        Sz = basis.quantize(gqcpy.ElectronicSpin_zOperator())
        self.Sz = Sz

        # S-S+.
        Smp_a = np.zeros_like(S_aa)
        Smp_b = S_bb
        S_mp = gqcpy.ScalarUSQOneElectronOperator_d(gqcpy.ScalarUSQOneElectronOperatorComponent_d(Smp_a), gqcpy.ScalarUSQOneElectronOperatorComponent_d(Smp_b))

        # S_z^2 operator.
        Sz2_a = 0.25 * S_aa
        Sz2_b = 0.25 * S_bb
        Sz2 = gqcpy.ScalarUSQOneElectronOperator_d(gqcpy.ScalarUSQOneElectronOperatorComponent_d(Sz2_a), gqcpy.ScalarUSQOneElectronOperatorComponent_d(-Sz2_b))

        # Total one-electron component.
        spin = Sz + S_mp + Sz2

        # Mixed two-electron component.
        N = int(basis.numberOfSpinors()/2)
        tensor_aa = gqcpy.ScalarPureUSQTwoElectronOperatorComponent_d(np.zeros((N, N, N, N)))
        tensor_ab = gqcpy.ScalarMixedUSQTwoElectronOperatorComponent_d(-1*np.einsum('pr,qs->pqrs', S_aa, S_bb))
        tensor_ba = gqcpy.ScalarMixedUSQTwoElectronOperatorComponent_d(-1*np.einsum('pr,qs->pqrs', S_bb, S_aa))
        tensor_bb = gqcpy.ScalarPureUSQTwoElectronOperatorComponent_d(np.zeros((N, N, N, N)))
        tensor = gqcpy.ScalarUSQTwoElectronOperator_d(tensor_aa, tensor_ab, tensor_ba, tensor_bb)

        self.matrix_operator = spin
        self.tensor_operator = tensor

        
        # For later uses, the operator object contains the "name" of the resulting expectation values.
        self.operator_result = "S^2 value"



class ConstrainedGHF(UnresolvedOptimizer):

    def __init__ (self, molecule, basis_set, operator, convergence_threshold=1.0e-06, DIIS=False):
        # Check compatibility of the operator type.
        assert (operator.type() == "UnresolvedSpinorBasisOperator"), "Only `UnresolvedSpinorBasisOperator` can be used with `ConstrainedGHF`."

        # Initialize the GHF spinor basis.
        spinor_basis = gqcpy.GSpinorBasis_d(molecule, basis_set)

        # This basis can now quantize the Hamiltonian.
        self.sq_hamiltonian = spinor_basis.quantize(gqcpy.FQMolecularHamiltonian(molecule))

        # We will need the number of electrons, as well as the total number of Spinors later on.        
        self.K = spinor_basis.numberOfSpinors()
        self.N = molecule.numberOfElectrons()

        # Save the overlap and nuclear repulsion operators. 
        self.nuclear_repulsion = gqcpy.NuclearRepulsionOperator(molecule.nuclearFramework()).value()
        self.overlap = spinor_basis.quantize(gqcpy.OverlapOperator())

        # Save the operator you want to constrain.
        self.operator = operator(spinor_basis)

        # Select the type of solver for the SCF algorithm.
        if DIIS:
            self.solver = "DIIS"
        else:
            self.solver = "Plain"

        self.convergence_threshold = convergence_threshold


    # A solver for the GHF problem.
    # The GHF solver always takes a random guess, in order to activate the `off-diagonal` blocks. 
    def _solveGHFProblem(self, hamiltonian):
        # To solve the GHF problem we need an environment and a solver.
        environment = gqcpy.GHFSCFEnvironment_d.WithCoreGuess(self.N, hamiltonian, self.overlap)
        
        if self.solver == "Plain":
            solver = gqcpy.GHFSCFSolver_d.Plain(threshold=self.convergence_threshold, maximum_number_of_iterations=250000)
        elif self.solver == "DIIS":
            solver = gqcpy.GHFSCFSolver_d.DIIS(minimum_subspace_dimension=6, maximum_subspace_dimension=6, threshold=self.convergence_threshold, maximum_number_of_iterations=250000)

        qc_structure = gqcpy.GHF_d.optimize(solver, environment)

        # For generalized Hartree-Fock, a stability check is always performed.
        # Transform the hamiltonian to MO basis and calculate the stability matrices. Print the resulting stabilities of the wavefunction model.
        # coefficients = qc_structure.groundStateParameters().expansion()
        # MO_hamiltonian = hamiltonian.transformed(coefficients)
        # stability_matrices = qc_structure.groundStateParameters().calculateStabilityMatrices(MO_hamiltonian)

        # internal_stability = stability_matrices.isInternallyStable(-1e-5)

        # while internal_stability is False:
        #     print("**************************************************************")
        #     print("There is an internal instability. Follow it using the Hessian.")
        #     print("**************************************************************")

        #     # Rotate the coefficients in the direction of the lowest Hessian eigenvector.
        #     rotation = stability_matrices.instabilityRotationMatrix(self.N, self.K-self.N)
        #     coefficients_rotated = coefficients.rotated(rotation)

        #     # Perform a new SCF calculation with the rotated coefficients as initial guess.
        #     environment_rotated = gqcpy.GHFSCFEnvironment_d(self.N, hamiltonian, self.overlap, coefficients_rotated)
        #     qc_structure = gqcpy.GHF_d.optimize(solver, environment_rotated)
        #     coefficients_2 = qc_structure.groundStateParameters().expansion()

        #     # Perform a new stability check. Print the resulting stabilities.
        #     hamiltonian_MO_2 = hamiltonian.transformed(coefficients_2)
        #     stability_matrices_2 = qc_structure.groundStateParameters().calculateStabilityMatrices(hamiltonian_MO_2)

        #     # Print the new stability consitions.
        #     stability_matrices_2.printStabilityDescription()

        #     # Update the internal stability parameter.
        #     internal_stability = stability_matrices_2.isInternallyStable(-1e-5)

        #     if internal_stability:
        #         print("**************************************************************")

        return qc_structure.groundStateEnergy(), qc_structure.groundStateParameters()


    # A function to calculate the energy and expectation value value of GHF at a certain multiplier. 
    def calculateEnergyAndExpectationValue(self, multiplier, return_parameters=False, verbose=0):
        modified_hamiltonian = self.sq_hamiltonian - (multiplier * self.operator.matrix_operator) - (multiplier * self.operator.tensor_operator)
        gs_energy, gs_parameters = self._solveGHFProblem(modified_hamiltonian)

        # Calculate the expectation value of the operator.
        OneDM = gs_parameters.calculateScalarBasis1DM()
        TwoDM = gs_parameters.calculateScalarBasis2DM()
        expectation_value = self.operator.matrix_operator.calculateExpectationValue(OneDM)[0] + self.operator.tensor_operator.calculateExpectationValue(TwoDM)[0]

        # Perform the energy correction
        energy = gs_energy + ((multiplier * expectation_value) + self.nuclear_repulsion)

        # Print the progress of which mu values have been completed if verbose >= 2.
        if verbose >= 2:
            print("--------------------------------------------------------")
            print("Mu = " + str(np.around(multiplier, 2)) + " done.")
            print("Sz = ", self.operator.Sz.calculateExpectationValue(OneDM)[0])

        if return_parameters:
            return energy, expectation_value, gs_parameters
        else:
            return energy, expectation_value
        

class ConstrainedUHF(UnresolvedOptimizer):

    def __init__ (self, molecule, number_of_alpha_electrons, number_of_beta_electrons, basis_set, operator, convergence_threshold=1.0e-08, DIIS=False, stability=False, random_guess=False):
        # Check compatibility of the operator type.
        assert (operator.type() == "ResolvedSpinorBasisOperator"), "Only `ResolvedSpinorBasisOperator` can be used with `ConstrainedUHF`."

        # Initialize the UHF spin orbital basis.
        spinor_basis = gqcpy.USpinOrbitalBasis_d(molecule, basis_set)

        # This basis can now quantize the Hamiltonian.
        self.sq_hamiltonian = spinor_basis.quantize(gqcpy.FQMolecularHamiltonian(molecule))

        # We will need the number of electrons, as well as the total number of Spinors later on.        
        self.Ka = spinor_basis.numberOfSpinors() // 2
        self.Kb = spinor_basis.numberOfSpinors() // 2

        self.Na = number_of_alpha_electrons
        self.Nb = number_of_beta_electrons

        # Save the overlap and nuclear repulsion operators. 
        self.nuclear_repulsion = gqcpy.NuclearRepulsionOperator(molecule.nuclearFramework()).value()
        self.overlap = spinor_basis.quantize(gqcpy.OverlapOperator())

        # Save the operator you want to constrain.
        self.operator = operator(spinor_basis)

        # Select the type of solver for the SCF algorithm.
        if DIIS:
            self.solver = "DIIS"
        else:
            self.solver = "Plain"

        self.stability = stability
        self.random_guess= random_guess

        self.convergence_threshold = convergence_threshold


    # A solver for the UHF problem.
    # The UHF solver always takes a random, different guess for alpha and beta, in order to break away from the RHF solution. 
    def _solveUHFProblem(self, hamiltonian):
        # Generate a random guess in order to find a true UHF solution.
        if self.random_guess:
            rand.seed(2)
            random_matrix_alpha = np.random.rand(self.Ka, self.Ka)
            random_matrix_alpha_transpose = random_matrix_alpha.T
            symmetric_random_matrix_alpha = random_matrix_alpha + random_matrix_alpha_transpose
            _, alpha_guess = np.linalg.eigh(symmetric_random_matrix_alpha)

            rand.seed(3)
            random_matrix_beta = np.random.rand(self.Kb, self.Kb)
            random_matrix_beta_transpose = random_matrix_beta.T
            symmetric_random_matrix_beta = random_matrix_beta + random_matrix_beta_transpose
            _, beta_guess = np.linalg.eigh(symmetric_random_matrix_beta)

            guess = gqcpy.UTransformation_d(gqcpy.UTransformationComponent_d(alpha_guess), gqcpy.UTransformationComponent_d(beta_guess))

            # To solve the GHF problem we need an environment and a solver.
            environment = gqcpy.UHFSCFEnvironment_d(self.Na, self.Nb, hamiltonian, self.overlap, guess)
        else:
            environment = gqcpy.UHFSCFEnvironment_d.WithCoreGuess(self.Na, self.Nb, hamiltonian, self.overlap)
        
        if self.solver == "Plain":
            solver = gqcpy.UHFSCFSolver_d.Plain(threshold=self.convergence_threshold, maximum_number_of_iterations=250000)
        elif self.solver == "DIIS":
            solver = gqcpy.UHFSCFSolver_d.DIIS(minimum_subspace_dimension=6, maximum_subspace_dimension=6, threshold=self.convergence_threshold, maximum_number_of_iterations=250000)

        qc_structure = gqcpy.UHF_d.optimize(solver, environment)

        # For unrestricted Hartree-Fock, a stability check is always performed.
        # Transform the hamiltonian to MO basis and calculate the stability matrices. Print the resulting stabilities of the wavefunction model.
        if self.stability:
            coefficients = qc_structure.groundStateParameters().expansion()
            MO_hamiltonian = hamiltonian.transformed(coefficients)
            stability_matrices = qc_structure.groundStateParameters().calculateStabilityMatrices(MO_hamiltonian)

            internal_stability = stability_matrices.isInternallyStable(-1e-5)

            while internal_stability is False:
                print("**************************************************************")
                print("There is an internal instability. Follow it using the Hessian.")
                print("**************************************************************")

                # Rotate the coefficients in the direction of the lowest Hessian eigenvector.
                rotation = stability_matrices.instabilityRotationMatrix(self.Na, self.Nb, self.Ka-self.Na, self.Kb-self.Nb)
                coefficients_rotated = coefficients.rotated(rotation)

                # Perform a new SCF calculation with the rotated coefficients as initial guess.
                environment_rotated = gqcpy.UHFSCFEnvironment_d(self.Na, self.Nb, hamiltonian, self.overlap, coefficients_rotated)
                qc_structure = gqcpy.UHF_d.optimize(solver, environment_rotated)
                coefficients_2 = qc_structure.groundStateParameters().expansion()

                # Perform a new stability check. Print the resulting stabilities.
                hamiltonian_MO_2 = hamiltonian.transformed(coefficients_2)
                stability_matrices_2 = qc_structure.groundStateParameters().calculateStabilityMatrices(hamiltonian_MO_2)

                # Print the new stability conditions.
                stability_matrices_2.printStabilityDescription()

                # Update the internal stability parameter.
                internal_stability = stability_matrices_2.isInternallyStable(-1e-5)

                if internal_stability:
                    print("**************************************************************")

        return qc_structure.groundStateEnergy(), qc_structure.groundStateParameters()


    # A function to calculate the energy and expectation value value of UHF at a certain multiplier. 
    def calculateEnergyAndExpectationValue(self, multiplier, return_parameters=False, verbose=0):
        # Modify the Hamiltonian with the given multiplier.
        # Unrestricted behaves difficultly so we extract all the wanted operator components.
        # One electron
        alpha_mat = gqcpy.ScalarUSQOneElectronOperatorComponent_d(copy.deepcopy(self.operator.matrix_operator.alpha.parameters()) * multiplier)
        beta_mat = gqcpy.ScalarUSQOneElectronOperatorComponent_d(copy.deepcopy(self.operator.matrix_operator.beta.parameters()) * multiplier)

        # two electron.
        ab_tens = gqcpy.ScalarMixedUSQTwoElectronOperatorComponent_d(copy.deepcopy(self.operator.tensor_operator.alphaAlpha().parameters()) * multiplier)
        aa_tens = gqcpy.ScalarPureUSQTwoElectronOperatorComponent_d(copy.deepcopy(self.operator.tensor_operator.alphaBeta().parameters()) * multiplier)
        bb_tens = gqcpy.ScalarPureUSQTwoElectronOperatorComponent_d(copy.deepcopy(self.operator.tensor_operator.betaAlpha().parameters()) * multiplier)
        ba_tens = gqcpy.ScalarMixedUSQTwoElectronOperatorComponent_d(copy.deepcopy(self.operator.tensor_operator.betaBeta().parameters()) * multiplier)

        unrestricted_one_e = gqcpy.ScalarUSQOneElectronOperator_d(alpha_mat, beta_mat)
        unrestricted_two_e = gqcpy.ScalarUSQTwoElectronOperator_d(aa_tens, ab_tens, ba_tens, bb_tens)

        # modify the hamiltonian.
        modified_hamiltonian = self.sq_hamiltonian - unrestricted_one_e - unrestricted_two_e

        # Run the HF calulation.
        gs_energy, gs_parameters = self._solveUHFProblem(modified_hamiltonian)

        # Calculate the expectation value of the operator.
        #print(gs_parameters.calculateScalarBasis1DM().alpha.matrix())
        OneDM = gs_parameters.calculateScalarBasis1DM()
        TwoDM = gs_parameters.calculateScalarBasis2DM()

        total_expectation_value = self.operator.matrix_operator.calculateExpectationValue(OneDM)[0] + self.operator.tensor_operator.calculateExpectationValue(TwoDM)[0]

        # Calculate the energy by correcting the ground state energy of the modified Hamiltonian.
        energy = gs_energy + (multiplier * total_expectation_value) + self.nuclear_repulsion

        # Print the progress of which mu values have been completed if verbose >= 2.
        if verbose >= 2:
            print("--------------------------------------------------------")
            print("Multiplier: " + str(np.around(multiplier, 2)) + " done.")
            print("Sz = ", np.around(self.operator.Sz.calculateExpectationValue(OneDM)[0], 2))
        if return_parameters:
            return energy, total_expectation_value, gs_parameters
        else:
            return energy, total_expectation_value

def CUHF_state(S_Squared, molecules, Na, Nb, basis_set, method, bracket=[-1, 1]):
    CUHF = []

    for mol in molecules:
        if S_Squared == 0.0:
            object = ConstrainedUHF(mol, Na, Nb, basis_set, RealSSquaredUOperator, 1e-3, random_guess=False)
        else:
            object = ConstrainedUHF(mol, Na, Nb, basis_set, RealSSquaredUOperator, 1e-3, random_guess=True)
        if method == "golden":
            output = GQCC.ExpectationValueSearch.unresolved(object, GQCC.UnresolvedOptimizationFunctions.GoldenLineSearch, [S_Squared], return_parameters=True, threshold=1e-3, verbose=2)
        elif method == "bisect":
            output = GQCC.ExpectationValueSearch.unresolved(object, GQCC.UnresolvedOptimizationFunctions.Bisect, [S_Squared], return_parameters=True, threshold=1e-3, verbose=2, bracket=bracket, tolerance=1e-3)
        output['alpha_parameters'] = output.apply(lambda row: getAlphaMatrix(row["parameters"]), axis=1)
        output['beta_parameters'] = output.apply(lambda row: getBetaMatrix(row["parameters"]), axis=1)
        CUHF.append(output)
    
    return CUHF