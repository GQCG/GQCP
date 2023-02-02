import numpy as np
import gqcpy

def UHF(molecule, N_a, N_b, basis_set, stability=True):

    spinor_basis = gqcpy.USpinOrbitalBasis_d(molecule, basis_set)

    fq_ham = gqcpy.FQMolecularHamiltonian(molecule)
    usq_ham = spinor_basis.quantize(fq_ham)
    S = spinor_basis.quantize(gqcpy.OverlapOperator())
    nuc_rep = gqcpy.NuclearRepulsionOperator(molecule.nuclearFramework()).value()    

    environment = gqcpy.UHFSCFEnvironment_d.WithCoreGuess(N_a, N_b, usq_ham, S)
    solver = gqcpy.UHFSCFSolver_d.Plain(threshold=1.0e-08, maximum_number_of_iterations=10000)   
    qc_structure = gqcpy.UHF_d.optimize(solver, environment)

    K_alpha = spinor_basis.numberOfSpinors() // 2
    K_beta = K_alpha

    Oa, Va = N_a, K_alpha - N_a
    Ob, Vb = N_b, K_beta - N_b

    coefficients = qc_structure.groundStateParameters().expansion()
    MO_hamiltonian = usq_ham.transformed(coefficients)
    stability_matrices = qc_structure.groundStateParameters().calculateStabilityMatrices(MO_hamiltonian)
    internal_stability = stability_matrices.isInternallyStable(-1e-5)

    if stability:
        while internal_stability is False:
            print("*********************")
            print("following instability")
            # Rotate the coefficients in the direction of the lowest Hessian eigenvector.
            rotation = stability_matrices.instabilityRotationMatrix(Oa, Ob, Va, Vb)
            coefficients_rotated = coefficients.rotated(rotation)

            # Perform a new SCF calculation with the rotated coefficients as initial guess.
            environment_rotated = gqcpy.UHFSCFEnvironment_d(N_a, N_b, usq_ham, S, coefficients_rotated)
            qc_structure = gqcpy.UHF_d.optimize(solver, environment_rotated)
            coefficients_2 = qc_structure.groundStateParameters().expansion()

            # Perform a new stability check. Print the resulting stabilities.
            hamiltonian_MO_2 = usq_ham.transformed(coefficients_2)
            stability_matrices_2 = qc_structure.groundStateParameters().calculateStabilityMatrices(hamiltonian_MO_2)
            stability_matrices_2.printStabilityDescription()
            internal_stability = stability_matrices_2.isInternallyStable(-1e-5)
        
    print("UHF energy: ", qc_structure.groundStateEnergy() + nuc_rep)
    
    return qc_structure.groundStateEnergy() + nuc_rep, qc_structure.groundStateParameters()


def UHFPoint(molecule, N_a, N_b, basis_set="STO-3G"):

    print("N_alpha, N_beta: {}, {}".format(N_a, N_b))
    spinor_basis = gqcpy.USpinOrbitalBasis_d(molecule, basis_set)
    fq_ham = gqcpy.FQMolecularHamiltonian(molecule)
    usq_ham = spinor_basis.quantize(fq_ham)

    K_alpha = spinor_basis.numberOfSpinors() // 2
    K_beta = K_alpha
    print("K_alpha, K_beta: {}, {}".format(K_alpha, K_beta))
    S = spinor_basis.quantize(gqcpy.OverlapOperator())
    nuc_rep = gqcpy.NuclearRepulsionOperator(molecule.nuclearFramework()).value()    

    environment = gqcpy.UHFSCFEnvironment_d.WithCoreGuess(N_a, N_b, usq_ham, S)
    solver = gqcpy.UHFSCFSolver_d.Plain(threshold=1.0e-07, maximum_number_of_iterations=250000)


    qc_structure = gqcpy.UHF_d.optimize(solver, environment)
    print(solver.numberOfIterations())
    Oa, Va = N_a, K_alpha - N_a
    Ob, Vb = N_b, K_beta - N_b
    print("(Oa, Va): ({},{})".format(Oa, Va))
    print("(Ob, Vb): ({},{})".format(Ob, Vb))


    # Stability check
    coefficients = qc_structure.groundStateParameters().expansion()
    MO_hamiltonian = usq_ham.transformed(coefficients)
    stability_matrices = qc_structure.groundStateParameters().calculateStabilityMatrices(MO_hamiltonian)
    internal_stability = stability_matrices.isInternallyStable(-1e-5)

    while internal_stability is False:
        print("*********************")
        print("following instability")
        # Rotate the coefficients in the direction of the lowest Hessian eigenvector.
        rotation = stability_matrices.instabilityRotationMatrix(Oa, Ob, Va, Vb)
        coefficients_rotated = coefficients.rotated(rotation)

        # Perform a new SCF calculation with the rotated coefficients as initial guess.
        environment_rotated = gqcpy.UHFSCFEnvironment_d(N_a, N_b, usq_ham, S, coefficients_rotated)
        qc_structure = gqcpy.UHF_d.optimize(solver, environment_rotated)
        coefficients_2 = qc_structure.groundStateParameters().expansion()

        # Perform a new stability check. Print the resulting stabilities.
        hamiltonian_MO_2 = usq_ham.transformed(coefficients_2)
        stability_matrices_2 = qc_structure.groundStateParameters().calculateStabilityMatrices(hamiltonian_MO_2)
        stability_matrices_2.printStabilityDescription()
        internal_stability = stability_matrices_2.isInternallyStable(-1e-5)
    print("    distance done    ")
    print(molecule)
    print("*********************")
    return qc_structure.groundStateEnergy() + nuc_rep

