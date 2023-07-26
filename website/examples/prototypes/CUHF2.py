import numpy as np
import scipy.linalg as la
import gqcpy

class CUHF2:

    def __init__(self, molecule, basis_set, Na, Nb):

        # Necessary variables
        spinor_basis = gqcpy.USpinOrbitalBasis_d(molecule, basis_set)
        usq_ham = spinor_basis.quantize(gqcpy.FQMolecularHamiltonian(molecule))


        # Integrals.
        self.core_alpha = usq_ham.core().alpha.parameters()
        self.core_beta = usq_ham.core().beta.parameters()

        self.eri_alpha = usq_ham.twoElectron().alphaAlpha().parameters()
        self.eri_beta = usq_ham.twoElectron().betaBeta().parameters()

        self.S = spinor_basis.quantize(gqcpy.OverlapOperator()).alpha.parameters()
        self.nuc_rep = gqcpy.NuclearRepulsionOperator(molecule.nuclearFramework()).value()

        # Information.
        self.dimension = np.shape(self.S)[0]
        self.Na = Na
        self.Nb = Nb

        # SCF environment variables.
        self.alpha_density = None
        self.beta_density = None

        self.Dplus = None
        self.Dmin = None

        self.alpha_fock_matrix = None
        self.beta_fock_matrix = None
        self.total_fock_matrix = None

        self.modified_fock_alpha = None
        self.modified_fock_beta = None

        self.alpha_coefficients = None
        self.beta_coefficients = None

        self.energy = 0
        self.converged = False
        self.iter = 0


    def _calculateGuessDensities(self):
        np.random.seed(2)
        xa = np.random.rand(self.dimension, self.dimension)
        xa_t = xa.T
        ya = xa + xa_t
        _, self.alpha_density = la.eigh(ya)
        #print(self.alpha_density)

        np.random.seed(3)
        xb = np.random.rand(self.dimension, self.dimension)
        xb_t = xb.T
        yb = xb + xb_t
        _, self.beta_density = la.eigh(yb)
        #print(self.beta_density)

    def _calculateDPlusMin(self):
        self.Dplus = self.alpha_density + self.beta_density
        self.Dmin = self.alpha_density - self.beta_density        

    def _calculateTotalFockMatrix(self):
        Jt = np.einsum('pqrs,rs->pq', self.eri_alpha, self.Dplus, optimize=True)
        Kt = np.einsum('prqs,rs->pq', self.eri_alpha, self.Dplus, optimize=True)
        self.total_fock_matrix = self.core_alpha + Jt - (0.5 * Kt)

    def _calculateAlphaBetaFockMatrices(self):
        Ks = np.einsum('prqs,rs->pq', self.eri_alpha, self.Dmin, optimize=True)
        self.alpha_fock_matrix = self.total_fock_matrix - (0.5 * Ks)
        self.beta_fock_matrix = self.total_fock_matrix + (0.5 * Ks)

    def _calculateModifiedFockMatrices(self, mu):
        self.modified_fock_alpha = self.alpha_fock_matrix - 2 * mu * (self.S @ self.beta_density @ self.S)
        self.modified_fock_beta = self.beta_fock_matrix - 2 * mu * (self.S @ self.alpha_density @ self.S)

    def _calculateAlphaBetaDensity(self):
        _, self.alpha_coefficients = la.eigh(self.modified_fock_alpha, self.S)
        occ_a = self.alpha_coefficients[:, :self.Na]
        self.alpha_density = np.einsum("ij, kj -> ik", occ_a, occ_a.conj())

        _, self.beta_coefficients = la.eigh(self.modified_fock_beta, self.S)
        occ_b = self.beta_coefficients[:, :self.Nb]
        self.beta_density = np.einsum("ij, kj -> ik", occ_b, occ_b.conj())

    def _calculateElectronicEnergy(self):
        # Za = self.core_alpha + self.modified_fock_alpha
        # Zb = self.core_beta + self.modified_fock_beta
        Za = self.core_alpha + self.alpha_fock_matrix
        Zb = self.core_beta + self.beta_fock_matrix

        E_alpha = 0.5 * np.einsum("ij, ij->", Za, self.alpha_density)
        E_beta = 0.5 * np.einsum("ij, ij->", Zb, self.beta_density)

        return E_alpha + E_beta

    def SCF(self, mu, maxiter=1000, threshold=1e-8):
        
        # Start from two different alpha and beta guesses.
        self._calculateGuessDensities()

        while self.converged is False and self.iter < maxiter:
            
            # Set up the Dplus and Dmin matrices
            self._calculateDPlusMin()

            # Calculate the modified fock matrices
            self._calculateTotalFockMatrix()
            self._calculateAlphaBetaFockMatrices()
            self._calculateModifiedFockMatrices(mu)

            # Calculate the energy.
            new_energy = self._calculateElectronicEnergy()

            # Calculate new densityu matrices.
            self._calculateAlphaBetaDensity()

            # Check for convergence and remain under the maximally allowed iterations.
            if np.abs(new_energy - self.energy) < threshold:
                self.iter += 1
                self.converged = True
                print("converged!")
            else:
                self.iter += 1
                self.energy = new_energy
            
            if self.iter == maxiter:
                print("not converged :(")
                self.energy = np.nan
                return

    def calculateKValue(self):
        k = np.trace(self.S @ self.beta_density @ self.S @ self.alpha_density)
        return k

    def calculateSSquared(self):
        Sz = (self.Na - self.Nb) / 2
        S2 = Sz * (Sz + 1) + self.Nb - self.calculateKValue()
        return S2







        

    
