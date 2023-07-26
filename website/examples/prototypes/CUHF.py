import numpy as np
import scipy.linalg as la
import gqcpy


def CUHF(molecule, N_a, N_b, basis_set, mu, log=False):

    spinor_basis = gqcpy.USpinOrbitalBasis_d(molecule, basis_set)  # in AO basis; atomic spin-orbitals
    S = spinor_basis.overlap()  # in AO basis

    sq_hamiltonian = spinor_basis.quantize(gqcpy.FQMolecularHamiltonian(molecule))  # in AO basis

    environment = gqcpy.UHFSCFEnvironment_d.WithCoreGuess(N_a, N_b, sq_hamiltonian, S)
    solver = gqcpy.UHFSCFSolver_d.Plain(threshold=1.0e-04, maximum_number_of_iterations=25000)

    def constrain_function(environment):
        def two_electron(da, db):
            j_a = np.einsum('pqrs,rs->pq', sq_hamiltonian.twoElectron().alphaAlpha().parameters(), da, optimize=True)
            j_b = np.einsum('pqrs,rs->pq', sq_hamiltonian.twoElectron().betaBeta().parameters(), db, optimize=True)
            k_a = np.einsum('prqs,rs->pq', sq_hamiltonian.twoElectron().alphaAlpha().parameters(), da, optimize=True)
            k_b = np.einsum('prqs,rs->pq', sq_hamiltonian.twoElectron().betaBeta().parameters(), db, optimize=True)
            return j_a, j_b, k_a, k_b

        # Ca = environment.coefficient_matrices[-1].alpha.matrix()
        Ca = spinor_basis.lowdinOrthonormalization().alpha.matrix()
        d_a = environment.density_matrices[-1].alpha.matrix()
        d_b = environment.density_matrices[-1].beta.matrix()
        j_a, j_b, k_a, k_b = two_electron(d_a, d_b)

        f_p = 0.5 * (2 * (j_a + j_b) - k_a - k_b)
        f_m = -0.5 * (k_a - k_b)

        p = (d_a + d_b) / 2
        p = la.inv(Ca) @ p @ la.inv(Ca.T)
        nat_occ_num, nat_occ_vec = la.eigh(p)
        nat_occ_vec = np.flip(nat_occ_vec, axis=1)

        f_m = nat_occ_vec.T @ Ca.T @ f_m @ Ca @ nat_occ_vec

        f_m_NO = np.copy(f_m)
        f_m_NO[:N_b, N_a:] = 0.0
        f_m_NO[N_a:, :N_b] = 0.0
        f_m_NO = np.linalg.inv(Ca.T) @ nat_occ_vec @ f_m_NO @ nat_occ_vec.T @ np.linalg.inv(Ca)
        
        f_a = sq_hamiltonian.core().alpha.parameters() + f_p + f_m_NO
        f_b = sq_hamiltonian.core().beta.parameters() + f_p - f_m_NO

        spin_contamination_matrix_a = f_a - environment.fock_matrices[-1].alpha.parameters()
        spin_contamination_matrix_b = environment.fock_matrices[-1].beta.parameters() - f_b

        fock_con_a = environment.fock_matrices[-1].alpha.parameters() + mu * spin_contamination_matrix_a
        fock_con_b = environment.fock_matrices[-1].beta.parameters() - mu * spin_contamination_matrix_b

        new_fock_matrices = gqcpy.ScalarUSQOneElectronOperator_d(gqcpy.ScalarUSQOneElectronOperatorComponent_d(fock_con_a), gqcpy.ScalarUSQOneElectronOperatorComponent_d(fock_con_b))
    
        environment.replace_current_fock_matrix(new_fock_matrices)
    # def constrain_function(environment):

    #     F = environment.fock_matrices[-1]
    #     D = environment.density_matrices[-1]
    #     N = environment.N  # Number of alpha and beta electrons.

    #     # Form the closed-shell Fock matrix and the UHF modification in the AO basis.
    #     F_cs = (F.alpha + F.beta) / 2.0
    #     Delta_UHF = (F.alpha - F.beta) / 2.0
    
    #     # Form the modified Fock matrices.
    #     F_aa = F_cs + Delta_UHF
    #     F_bb = F_cs - Delta_UHF

    #     # Find the natural occupation numbers and vectors by diagonalizing the charge-density matrix
    #     # in an orthonormal basis. For the orthonormal basis, we use the Löwdin basis of the alpha basis functions.
    #     X = spinor_basis.lowdinOrthonormalization().alpha  # X transforms from the AO basis to the Löwdin basis.
    
    #     P_AO = (D.alpha + D.beta) / 2.0
    #     P_MO = P_AO.transformed(X)

    #     natural_occupation_numbers, V = np.linalg.eigh(P_MO.matrix())  # Use NumPy for the diagonalization.

    #     # Construct the Langrange multipliers to add them to the 'constrained' Fock matrices. They
    #     # should be constructed in the basis of the natural occupations.
    #     V = gqcpy.UTransformationComponent_d(V)  # V transforms from the Löwdin basis to the natural occupations.

    #     Delta_UHF_NO = Delta_UHF.transformed(X).transformed(V)  # In the natural occupation (NO) basis.
    #     Delta_UHF_NO = Delta_UHF_NO.parameters()

    #     Lambda_NO = np.zeros(np.shape(Delta_UHF_NO))
    #     Lambda_NO[:N_b, N_a:] = -Delta_UHF_NO[:N_b, N_a:]
    #     Lambda_NO[N_a:, :N_b] = -Delta_UHF_NO[N_a:, :N_b]
    #     Lambda_NO = gqcpy.ScalarUSQOneElectronOperatorComponent_d(Lambda_NO)
    #     Lambda_AO = Lambda_NO.transformed(V.inverse()).transformed(X.inverse())  # In the AO basis.
    
    #     # Overwrite the most recent UHF Fock matrices with the CUHF modifications.
    #     F_alpha_constrained = F_aa + Lambda_AO
    #     F_beta_constrained = F_bb - Lambda_AO

    #     new_fock_matrices = gqcpy.ScalarUSQOneElectronOperator_d(F_alpha_constrained, F_beta_constrained)
    
    #     environment.replace_current_fock_matrix(new_fock_matrices)

    constrain_step = gqcpy.FunctionalStep_UHFSCFEnvironment_d(constrain_function, description="Replace the alpha- and beta- UHF Fock matrices by their constrained counterparts, as explained in the CUHF 'thesis' algorithm.")
    solver.insert(constrain_step, 2)

    # def print_func(environment):
    #     print("electronic E: ", environment.electronic_energies[-1])
    # print_step = constrain_step = gqcpy.FunctionalStep_UHFSCFEnvironment_d(print_func, description="print E")

    # solver.insert(print_step, 5)

    cuhf_qc_structure = gqcpy.UHF_d.optimize(solver, environment)
    cuhf_energy = cuhf_qc_structure.groundStateEnergy()
    cuhf_parameters = cuhf_qc_structure.groundStateParameters()

    nuc_rep = gqcpy.NuclearRepulsionOperator(molecule.nuclearFramework()).value()    
    print("CUHF energy: ", cuhf_energy + nuc_rep)

    return cuhf_energy + nuc_rep, cuhf_parameters
