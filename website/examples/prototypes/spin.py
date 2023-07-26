import numpy as np
import gqcpy

def SpinExpectationValues(molecule, Na, Nb, basis_set, parameters):
    occ_a = Na
    occ_b = Nb

    spinor_basis = gqcpy.USpinOrbitalBasis_d(molecule, basis_set)  # in AO basis; atomic spin-orbitals
    overlap = spinor_basis.overlap().alpha.parameters()  # in AO basis

    occ_indx_a = np.arange(occ_a)  # indices of the occupied alpha orbitals
    occ_indx_b = np.arange(occ_b)  # indices of the occupied beta orbitals

    occ_a_orb = parameters.expansion().alpha.matrix()[:, occ_indx_a]  # orbital coefficients associated with occupied alpha orbitals
    occ_b_orb = parameters.expansion().beta.matrix()[:, occ_indx_b]  # orbital coefficients associated with occupied beta orbitals

    s = occ_a_orb.conj().T @ overlap @ occ_b_orb  # Basically (alpha orbitals).T * S * (beta orbitals)

    ss_xy = (occ_a + occ_b) * 0.5 - np.einsum('ij,ij->', s.conj(), s)  # = S^2_x + S^2_y
    ss_z = (occ_b - occ_a)**2 * 0.25  # = S^2_z
    ss = (ss_xy + ss_z).real  # = S^2_total
    s_z = (occ_a - occ_b) / 2  # = S_z

    return s_z, ss

