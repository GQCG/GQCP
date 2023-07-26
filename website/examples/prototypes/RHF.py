import numpy as np
import gqcpy

def RHF(molecule, basis_set):

    spinor_basis = gqcpy.RSpinOrbitalBasis_d(molecule, basis_set)

    fq_ham = gqcpy.FQMolecularHamiltonian(molecule)
    sq_ham = spinor_basis.quantize(fq_ham)
    S = spinor_basis.quantize(gqcpy.OverlapOperator())
    nuc_rep = gqcpy.NuclearRepulsionOperator(molecule.nuclearFramework()).value()    

    environment = gqcpy.RHFSCFEnvironment_d.WithCoreGuess(molecule.numberOfElectrons(), sq_ham, S)
    solver = gqcpy.RHFSCFSolver_d.Plain(threshold=1.0e-08, maximum_number_of_iterations=25000)   
    rhf_objective = gqcpy.DiagonalRHFFockMatrixObjective_d(sq_ham)



    qc_structure = gqcpy.RHF_d.optimize(rhf_objective, solver, environment)
    print("RHF energy: ", qc_structure.groundStateEnergy() + nuc_rep)
    
    return qc_structure.groundStateEnergy() + nuc_rep, qc_structure.groundStateParameters()

def RHFPES(molecules, basis_set):
    RHF_E = []
    RHF_parameters = []
    RHF_parameters_resolved = []

    for mol in molecules:
        E, RHF_par = RHF(mol, basis_set)
        RHF_E.append(E)

    for parameters in RHF_parameters:
        expansion = parameters.expansion()
        resolved_expansion = gqcpy.UTransformation_d.FromRestricted(expansion)
        RHF_parameters_resolved.append(resolved_expansion)

    return RHF_E, RHF_parameters_resolved
