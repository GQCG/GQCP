import gqcpy
import numpy as np
import scipy.linalg as la

class UNOCI:

    def __init__(self, molecule, alpha, beta, basis_set, basis_state_vector):
        self.molecule = molecule
        self.Na = alpha
        self.Nb = beta
        self.spinor_basis = gqcpy.USpinOrbitalBasis_d(self.molecule, basis_set)
        fq_hamiltonian = gqcpy.FQMolecularHamiltonian(self.molecule)
        self.sq_hamiltonian = self.spinor_basis.quantize(fq_hamiltonian)
        self.S = self.spinor_basis.quantize(gqcpy.OverlapOperator())
        self.basis = gqcpy.UNonOrthogonalStateBasis_d(basis_state_vector, self.S, self.Na, self.Nb)


    def optimize(self):
        environment = gqcpy.NOCIEnvironment.Dense_d(self.sq_hamiltonian, self.basis, self.molecule)
        solver = gqcpy.GeneralizedEigenproblemSolver.Dense_d()
        qc_structure = gqcpy.NOCI_d(self.basis).optimize(solver, environment)
        return qc_structure.groundStateEnergy(), qc_structure.groundStateParameters()


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
        parameters = gqcpy.py_NOCIExpansion_UNonOrthogonalStateBasis_d(self.basis, gs_coefficients)

        return NOCI_val[0], parameters
