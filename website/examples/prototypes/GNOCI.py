import gqcpy
import scipy.linalg as la
import numpy as np

class GNOCI:

    def __init__(self, basis_state_vector, molecule, basis_set):
        self.molecule = molecule
        N = self.molecule.numberOfElectrons()
        self.spinor_basis = gqcpy.GSpinorBasis_d(self.molecule, basis_set)
        fq_hamiltonian = gqcpy.FQMolecularHamiltonian(self.molecule)
        self.sq_hamiltonian = self.spinor_basis.quantize(fq_hamiltonian)
        self.S = self.spinor_basis.quantize(gqcpy.OverlapOperator())
        self.basis = gqcpy.GNonOrthogonalStateBasis_d(basis_state_vector, self.S, N)


    def optimize(self):
        environment = gqcpy.NOCIEnvironment.Dense_d(self.sq_hamiltonian, self.basis, self.molecule)
        solver = gqcpy.GeneralizedEigenproblemSolver.Dense_d()
        qc_structure = gqcpy.NOCI_d(self.basis).optimize(solver, environment)
        return qc_structure.groundStateEnergy(), qc_structure.groundStateParameters()


    def optimize_excited_states(self, number_of_states):
        environment = gqcpy.NOCIEnvironment.Dense_d(self.sq_hamiltonian, self.basis, self.molecule)
        solver = gqcpy.GeneralizedEigenproblemSolver.Dense_d()
        qc_structure = gqcpy.NOCI_d(self.basis, number_of_states).optimize(solver, environment)

        energies = []
        parameters = []

        for i in range(number_of_states):
            energies.append(qc_structure.energy(i))
            parameters.append(qc_structure.parameters(i))

        return energies, parameters


    def optimize_full_C_output(self, number_of_basis_states):
        environment = gqcpy.NOCIEnvironment.Dense_d(self.sq_hamiltonian, self.basis, self.molecule)
        solver = gqcpy.GeneralizedEigenproblemSolver.Dense_d()
        qc_structure = gqcpy.NOCI_d(self.basis, number_of_states=number_of_basis_states).optimize(solver, environment)

        C_matrix = np.zeros((number_of_basis_states, number_of_basis_states))
        for i in range(number_of_basis_states):
            parameters = qc_structure.parameters(i)
            coeff = parameters.coefficients()
            C_matrix[:, [i]] = coeff.reshape((number_of_basis_states, 1))

        overlap = self.basis.evaluateOverlapOperator()

        return C_matrix, overlap


    def optimizeInSubspace(self, threshold):
        NOCI_Hamiltonian = self.basis.evaluateHamiltonianOperator(self.sq_hamiltonian, gqcpy.NuclearRepulsionOperator(self.molecule.nuclearFramework()))
        NOCI_overlap = self.basis.evaluateOverlapOperator()

        val, vec = la.eigh(NOCI_overlap)
        zero_index = 0

        for value in val:
            if value < threshold:
                zero_index += 1
        
        if zero_index > 0:
            print("The overlap eigenvalues are: ")
            print(val)
            print("Number of zero overlap values:")
            print(zero_index)

        vec = vec[:, zero_index:]
        subspace_Hamiltonian = vec.conjugate().T @ NOCI_Hamiltonian @ vec
        subspace_overlap = vec.conjugate().T @ NOCI_overlap @ vec

        if zero_index > 0:
            print("The new Hamiltonian dimension is: ")
            print(np.shape(subspace_Hamiltonian))
        
        # Solve the NOCI generalized eigenvalue problem.
        NOCI_val, NOCI_vec = la.eigh(subspace_Hamiltonian, subspace_overlap)

        # Transform the coefficients back to the original space.
        coefficients_backtransformed = vec @ NOCI_vec
        if zero_index > 0:
            print("The new NOCI coefficient matrix dimension is: ")
            print(np.shape(coefficients_backtransformed))

        gs_coefficients = coefficients_backtransformed[:, 0]
        parameters = gqcpy.py_NOCIExpansion_GNonOrthogonalStateBasis_d(self.basis, gs_coefficients)

        return NOCI_val[0], parameters


def NOCISz(ground_state_parameters, spinor_basis):
    Sz_operator = spinor_basis.quantize(gqcpy.ElectronicSpin_zOperator())
    DM = ground_state_parameters.calculate1DM()
    expectation_value = Sz_operator.calculateExpectationValue(DM)[0]
    return expectation_value
