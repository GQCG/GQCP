import gqcpy
import numpy as np

A = gqcpy.AdjacencyMatrix.Linear(4)

t = 1.0
U = 3.5
H = gqcpy.HoppingMatrix(A, t, U)

hubbard_hamiltonian = gqcpy.HubbardHamiltonian(H)

K = 4  # number of sites
N_P = 2  # number of electron pairs

onv_basis = gqcpy.SpinResolvedONVBasis(K, N_P, N_P)
solver = gqcpy.EigenproblemSolver.Dense()
environment = gqcpy.CIEnvironment.Dense(hubbard_hamiltonian, onv_basis)

energy = gqcpy.CI(onv_basis).optimize(solver, environment).groundStateEnergy()

np.testing.assert_almost_equal(energy, -2.135871608231944)
