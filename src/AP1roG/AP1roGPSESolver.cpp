#include "AP1roGPSESolver.hpp"


namespace GQCG {


/*
 * CONSTRUCTORS
 */
/**
 *  Constructor based on a given @param molecule and Hamiltonian parameters @param ham_par
 *
 *  The initial guess for the geminal coefficients is zero
 */
AP1roGPSESolver::AP1roGPSESolver(const GQCG::Molecule& molecule, const GQCG::HamiltonianParameters& ham_par) :
    K (ham_par.K),
    ham_par (ham_par),
    N_P (molecule.N / 2),
    initial_geminal_coefficients (GQCG::AP1roGGeminalCoefficients(this->N_P, this->K))
{
    // Check if we have an even number of electrons
    if ((molecule.N % 2) != 0) {
        throw std::invalid_argument("The given molecule has an odd number of electrons.");
    }
}



/*
 *  PUBLIC METHODS
 */
/**
 *  For a geminal coefficient g_mu, return its major index in the matrix of geminal coefficients.
 *
 *      Note that:
 *          - the major index is i (i.e. the subscript), since changes in i are not contiguous
 *          - i is in [0 ... N_P[
 */
size_t AP1roGPSESolver::matrixIndexMajor(size_t vector_index) const {

    return vector_index / (this->K - this->N_P);  // in the mathematical notes, we use the floor function, which is the same as integer division
}


/**
 *  For a geminal coefficient g_mu, return its minor index in the matrix of geminal coefficients.
 *
 *      Note that:
 *          - the minor index is a (i.e. the superscript), since changes in a are contiguous
 *          - a is in [N_P ... K[
 */
size_t AP1roGPSESolver::matrixIndexMinor(size_t vector_index) const {

    return vector_index % (this->K - this->N_P) + this->N_P;  // we add N_P since we want a to be in [N_P ... K[
}


/**
 *  For a geminal coefficient G_i^a, return its index in the vector of geminal coefficients.
 *
 *      Note that
 *          - i is in [0 ... N_P[       is the 'major' index (i.e. changes in i are not contiguous)
 *          - a is in [N_P ... K[       is the 'minor' index (i.e. changes in a are contiguous)
 */
size_t AP1roGPSESolver::vectorIndex(size_t i, size_t a) const {

    // Check for invalid values for i and a
    if (i >= this->N_P) {
        throw std::invalid_argument("i must be smaller than N_P.");
    }
    if (a < this->N_P) {
        throw std::invalid_argument("a must be larger than or equal to N_P.");
    }


    // The conversion from i and a to a single vector index is just a little more complicated than row-major storage.
    // If we were to use the row-major storage formula, we would end up with
    //      mu = a + (K - N_P) * i
    //
    // but since we would really like our indices abcd (virtuals) to start at N_P, we should subtract N_P accordingly
    return (a - this->N_P) + (this->K - this->N_P) * i;
}

}  // namespace GQCG
