import gqcpy
import numpy as np

def FCI(molecule, Na, Nb, basis_set, n=1):
    
    # Pre-requisites for a FCI calculation.
    N_a = Na
    N_b = Nb
        
    nuc_rep = gqcpy.NuclearRepulsionOperator(molecule.nuclearFramework()).value() # Calculate the nuclear repulsion value.
    spin_orbital_basis = gqcpy.RSpinOrbitalBasis_d(molecule, basis_set) # Generate the spin-orbital basis for the molecule in the given basis set.
    spin_orbital_basis.lowdinOrthonormalize() # For an FCI calculation, the spin orbital basis has to be orthonormal.
        
    sq_hamiltonian = spin_orbital_basis.quantize(gqcpy.FQMolecularHamiltonian(molecule)) # Create the second-quantized Hamiltonian of the molecule in the correct basis.
    K = spin_orbital_basis.numberOfSpatialOrbitals() # Determine the number of orbitals.
    S = spin_orbital_basis.quantize(gqcpy.OverlapOperator()) # Calculate the overlap matrix.
        
    onv_basis = gqcpy.SpinResolvedONVBasis(K, N_a, N_b) # Create the ONV basis of the correct size.
        
    # perform a dense FCI calculation.
    solver = gqcpy.EigenproblemSolver.Dense_d() # Choose the way you want to solve the FCI equations.
    environment = gqcpy.CIEnvironment.Dense(sq_hamiltonian, onv_basis) # Set up an environment in which you will perform the calculation.
    qc_structure = gqcpy.CI(onv_basis, number_of_states=n).optimize(solver, environment) 
                
    print("FCI energy: ", qc_structure.groundStateEnergy() + nuc_rep)
    
    if n == 1:
        return qc_structure.groundStateEnergy() + nuc_rep, qc_structure.groundStateParameters()
    if n > 1:
        expansions = []
        energies = []
        for i in range(n):
            # print("Excited state " + str(i) + ":", qc_structure.energy(i) + nuclear_repulsion)
            energies.append(qc_structure.energy(i) + nuc_rep)
            expansions.append(qc_structure.parameters(i))
        return energies, expansions
    

def FCIPES(molecules, Na, Nb, basis_set):
    FCI_energies = []
    for mol in molecules:
        FCI_E, _ = FCI(mol, Na, Nb, basis_set)
        FCI_energies.append(FCI_E)
    return FCI_energies
