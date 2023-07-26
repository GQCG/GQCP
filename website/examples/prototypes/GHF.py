import numpy as np
import numpy.random as rand
import gqcpy


def GHF(molecule, basis_set):

    spinor_basis = gqcpy.GSpinorBasis_d(molecule, basis_set)
    K = spinor_basis.numberOfSpinors()
    
    fq_ham = gqcpy.FQMolecularHamiltonian(molecule)
    sq_ham = spinor_basis.quantize(fq_ham)
    S = spinor_basis.quantize(gqcpy.OverlapOperator())
    nuc_rep = gqcpy.NuclearRepulsionOperator(molecule.nuclearFramework()).value()  

    # rand.seed(2)
    # random_matrix = np.random.rand(K, K)
    # random_matrix_transpose = random_matrix.T
    # symmetric_random_matrix = random_matrix + random_matrix_transpose
    # _, guess = np.linalg.eigh(symmetric_random_matrix)
  

    environment = gqcpy.GHFSCFEnvironment_d.WithCoreGuess(molecule.numberOfElectrons(), sq_ham, S)

    solver = gqcpy.GHFSCFSolver_d.Plain(threshold=1.0e-06, maximum_number_of_iterations=10000)   

    qc_structure = gqcpy.GHF_d.optimize(solver, environment)
    print("GHF energy: ", qc_structure.groundStateEnergy() + nuc_rep)
    
    return qc_structure.groundStateEnergy() + nuc_rep, qc_structure.groundStateParameters()