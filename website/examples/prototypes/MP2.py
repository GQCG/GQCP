import gqcpy
import numpy as np
import scipy.linalg as la
from .transform import *

class UMP2:

    def __init__(self, uhf_energy, uhf_parameters, hamiltonian):

        # The UHF data.
        self.UHF_energy = uhf_energy
        self.uhf_parameters = uhf_parameters
        self.hamiltonian = hamiltonian

        # The alpha and beta orbital energies.
        self.occupied_alpha_orbital_energies = np.array(uhf_parameters.occupiedOrbitalEnergies().alpha)
        self.occupied_beta_orbital_energies = np.array(uhf_parameters.occupiedOrbitalEnergies().beta)

        self.virtual_alpha_orbital_energies = np.array(uhf_parameters.virtualOrbitalEnergies().alpha)
        self.virtual_beta_orbital_energies = np.array(uhf_parameters.virtualOrbitalEnergies().beta)

        # The alpha and beta orbital coefficients, split in occupied/virtual.
        self.occupied_alpha_expansion = uhf_parameters.expansion().alpha.matrix()[:, :uhf_parameters.numberOfElectrons().alpha]
        self.occupied_beta_expansion = uhf_parameters.expansion().beta.matrix()[:, :uhf_parameters.numberOfElectrons().beta]

        self.virtual_alpha_expansion = uhf_parameters.expansion().alpha.matrix()[:, uhf_parameters.numberOfElectrons().alpha:]
        self.virtual_beta_expansion = uhf_parameters.expansion().beta.matrix()[:, uhf_parameters.numberOfElectrons().beta:]

        # Get the Hamiltonian properties.
        eri = hamiltonian.twoElectron().alphaAlpha().parameters()
        core_alpha = hamiltonian.core().alpha.parameters()
        core_beta = hamiltonian.core().beta.parameters()

        # Transform the two electron integrals.
        self.alpha_two_electron_integrals = tensorToMOBasis(eri, self.occupied_alpha_expansion, self.virtual_alpha_expansion, self.occupied_alpha_expansion, self.virtual_alpha_expansion).transpose(0, 2, 1, 3)
        self.beta_two_electron_integrals = tensorToMOBasis(eri, self.occupied_beta_expansion, self.virtual_beta_expansion, self.occupied_beta_expansion, self.virtual_beta_expansion).transpose(0, 2, 1, 3)
        self.mixed_two_electron_integrals = tensorToMOBasis(eri, self.occupied_alpha_expansion, self.virtual_alpha_expansion, self.occupied_beta_expansion, self.virtual_beta_expansion).transpose(0, 2, 1, 3)

        # Get the occupied alpha and beta Fock matrix parts.
        Ja = np.einsum('pqrs,rs->pq', eri, uhf_parameters.calculateScalarBasis1DM().alpha.matrix(), optimize=True)
        Ka = np.einsum('prqs,rs->pq', eri, uhf_parameters.calculateScalarBasis1DM().alpha.matrix(), optimize=True)
        Jb = np.einsum('pqrs,rs->pq', eri, uhf_parameters.calculateScalarBasis1DM().beta.matrix(), optimize=True)
        Kb = np.einsum('prqs,rs->pq', eri, uhf_parameters.calculateScalarBasis1DM().beta.matrix(), optimize=True)

        Fa = core_alpha + Ja - Ka
        Fb = core_beta + Jb - Kb

        self.occupied_alpha_fock_matrix = Fa[:uhf_parameters.numberOfElectrons().alpha, uhf_parameters.numberOfElectrons().alpha:]
        self.occupîed_beta_fock_matrix = Fb[:uhf_parameters.numberOfElectrons().beta, uhf_parameters.numberOfElectrons().beta:]

        # Calculate the MP2 denominators.
        e_denom_aa = self.occupied_alpha_orbital_energies.reshape(-1, 1, 1, 1)
        e_denom_aa = e_denom_aa - self.virtual_alpha_orbital_energies.reshape(-1, 1, 1)
        e_denom_aa = e_denom_aa + self.occupied_alpha_orbital_energies.reshape(-1, 1)
        e_denom_aa = e_denom_aa - self.virtual_alpha_orbital_energies
        self.e_denom_aa = 1 / e_denom_aa.transpose(0, 2, 1, 3)

        e_denom_bb = self.occupied_beta_orbital_energies.reshape(-1, 1, 1, 1)
        e_denom_bb = e_denom_bb - self.virtual_beta_orbital_energies.reshape(-1, 1, 1)
        e_denom_bb = e_denom_bb + self.occupied_beta_orbital_energies.reshape(-1, 1)
        e_denom_bb = e_denom_bb - self.virtual_beta_orbital_energies
        self.e_denom_bb = 1 / e_denom_bb.transpose(0, 2, 1, 3)

        e_denom_ab = self.occupied_alpha_orbital_energies.reshape(-1, 1, 1, 1)
        e_denom_ab = e_denom_ab - self.virtual_alpha_orbital_energies.reshape(-1, 1, 1)
        e_denom_ab = e_denom_ab + self.occupied_beta_orbital_energies.reshape(-1, 1)   
        e_denom_ab = e_denom_ab - self.virtual_beta_orbital_energies
        self.e_denom_ab = 1 / e_denom_ab.transpose(0, 2, 1, 3)


    def calculateMP2EnergyCorrection(self):
        # Calculate the alpha, same-spin contributions.
        E_ss_aa = np.einsum("ijab, ijab, ijab->", self.alpha_two_electron_integrals, self.alpha_two_electron_integrals - self.alpha_two_electron_integrals.swapaxes(2, 3), self.e_denom_aa, optimize=True)

        # calculate the beta, same-spin contributions.
        E_ss_bb = np.einsum("ijab, ijab, ijab->", self.beta_two_electron_integrals, self.beta_two_electron_integrals - self.beta_two_electron_integrals.swapaxes(2, 3), self.e_denom_bb, optimize=True)

        # calculate the opposite spin contributions.
        E_os = np.einsum("ijab, ijab, ijab->", self.mixed_two_electron_integrals, self.mixed_two_electron_integrals, self.e_denom_ab, optimize=True)

        # calculate and return the energy correction + the total MP2 energy.
        return (E_ss_aa + E_ss_bb + E_os), (self.UHF_energy + E_ss_aa + E_ss_bb + E_os)

    
    def _calculateT2Amplitudes(self, eri, denom, mixed=False):
        if mixed:
            return eri * denom
        else:
            anti_symmetrized = eri - eri.swapaxes(2, 3)
            return anti_symmetrized * denom


    def calculateMP2DensityMatrixCorrection(self):
        # Calculate T2 amplitudes.
        t2_aa = self._calculateT2Amplitudes(self.alpha_two_electron_integrals, self.e_denom_aa)
        t2_bb = self._calculateT2Amplitudes(self.beta_two_electron_integrals, self.e_denom_bb)
        t2_ab = self._calculateT2Amplitudes(self.mixed_two_electron_integrals, self.e_denom_ab, mixed=True)

        # calculate occupied and virtual MP2 densities.
        G_occupied_alpha = np.einsum('prqs,wrqs->pw', t2_aa.conj(), t2_aa) * -0.5
        G_occupied_alpha -= np.einsum('prqs,wrqs->pw', t2_ab.conj(), t2_ab)

        G_virtual_alpha = np.einsum('prqs,prvs->vq', t2_aa.conj(), t2_aa) * 0.5
        G_virtual_alpha += np.einsum('prqs,prvs->vq', t2_ab.conj(), t2_ab)

        G_occupied_beta = np.einsum('prqs,wrqs->pw', t2_bb.conj(), t2_bb) * -0.5
        G_occupied_beta -= np.einsum('rpqs,rwqs->pw', t2_ab.conj(), t2_ab)

        G_virtual_beta = np.einsum('prqs,prvs->vq', t2_bb.conj(), t2_bb) * 0.5
        G_virtual_beta += np.einsum('prsq,prsv->vq', t2_ab.conj(), t2_ab)

        # calculate the full MP2 alpha and beta density.
        # alpha part.
        D_occupied_alpha = G_occupied_alpha + G_occupied_alpha.conj().T
        D_virtual_alpha = G_virtual_alpha + G_virtual_alpha.conj().T

        D_alpha = np.zeros((self.uhf_parameters.numberOfSpinOrbitals(gqcpy.alpha), self.uhf_parameters.numberOfSpinOrbitals(gqcpy.alpha)))
        D_alpha[:self.uhf_parameters.numberOfElectrons().alpha, :self.uhf_parameters.numberOfElectrons().alpha] = D_occupied_alpha
        D_alpha[self.uhf_parameters.numberOfElectrons().alpha:, self.uhf_parameters.numberOfElectrons().alpha:] = D_virtual_alpha

        D_alpha *= 0.5
        D_alpha[np.diag_indices(self.uhf_parameters.numberOfElectrons().alpha)] += 1

        # Beta part.
        D_occupied_beta = G_occupied_beta + G_occupied_beta.conj().T
        D_virtual_beta = G_virtual_beta + G_virtual_beta.conj().T

        D_beta = np.zeros((self.uhf_parameters.numberOfSpinOrbitals(gqcpy.beta), self.uhf_parameters.numberOfSpinOrbitals(gqcpy.beta)))
        D_beta[:self.uhf_parameters.numberOfElectrons().beta, :self.uhf_parameters.numberOfElectrons().beta] = D_occupied_beta
        D_beta[self.uhf_parameters.numberOfElectrons().beta:, self.uhf_parameters.numberOfElectrons().beta:] = D_virtual_beta

        D_beta *= 0.5
        D_beta[np.diag_indices(self.uhf_parameters.numberOfElectrons().beta)] += 1

        eigvals_a, _ = la.eigh(D_alpha)
        eigvals_b, _ = la.eigh(D_beta)      

        print("N_alpha = ", np.sum(eigvals_a))
        print("N_beta = ", np.sum(eigvals_b))
        
        return D_alpha, D_beta

    
    # def calculateMP2CoefficientMatrixCorrection(self):
    #     # Calculate the density matrices.
    #     Da, Db = self.calculateMP2DensityMatrixCorrection()

    #     # Set up the Hamiltonian parameters.
    #     # Get the Hamiltonian properties.
    #     eri = self.hamiltonian.twoElectron().alphaAlpha().parameters()
    #     core_alpha = self.hamiltonian.core().alpha.parameters()
    #     core_beta = self.hamiltonian.core().beta.parameters()

    #     # Calculate new Fock matrices.
    #     # Get the occupied alpha and beta Fock matrix parts.
    #     Ja = np.einsum('pqrs,rs->pq', eri, Da, optimize=True)
    #     Ka = np.einsum('prqs,rs->pq', eri, Da, optimize=True)
    #     Jb = np.einsum('pqrs,rs->pq', eri, Db, optimize=True)
    #     Kb = np.einsum('prqs,rs->pq', eri, Db, optimize=True)

    #     Fa = core_alpha + Ja - Ka
    #     Fb = core_beta + Jb - Kb

    #     # Diagonalize the fock matrices to retrieve the coefficients.
    #     _, Ca = la.eigh(Fa)
    #     _, Cb = la.eigh(Fb)

    #     return Ca, Cb
