#ifndef AP1roGPSESolver_hpp
#define AP1roGPSESolver_hpp


#include "HamiltonianParameters/HamiltonianParameters.hpp"
#include "Molecule.hpp"
#include "AP1roG/AP1roGGeminalCoefficients.hpp"
#include "AP1roG/AP1roG.hpp"

namespace GQCG {

/**
 *
 */
class AP1roGPSESolver {
private:
    const size_t K;  // the number of special orbitals
    const size_t N_P;  // the number of electron pairs
    const GQCG::AP1roGGeminalCoefficients initial_geminal_coefficients;
    
    GQCG::HamiltonianParameters ham_par;

    GQCG::AP1roG solution;


public:
    // CONSTRUCTORS
    /**
     *  Constructor based on a given @param molecule and Hamiltonian parameters @param ham_par
     *
     *  The initial guess for the geminal coefficients is zero
     */
    AP1roGPSESolver(const GQCG::Molecule& molecule, const GQCG::HamiltonianParameters& ham_par);


    // GETTERS
    GQCG::AP1roG get_solution() const { return this->solution; }


    // PUBLIC METHODS
    /**
     *  For a geminal coefficient g_mu, return its major index in the matrix of geminal coefficients.
     *
     *      Note that:
     *          - the major index is i (i.e. the subscript), since changes in i are not contiguous
     *          - i is in [0 ... N_P[
     */
    size_t matrixIndexMajor(size_t vector_index) const;

    /**
     *  For a geminal coefficient g_mu, return its minor index in the matrix of geminal coefficients.
     *
     *      Note that:
     *          - the minor index is a (i.e. the superscript), since changes in a are contiguous
     *          - a is in [N_P ... K[
     */
    size_t matrixIndexMinor(size_t vector_index) const;

    /**
     *  For a geminal coefficient G_i^a, return its index in the vector of geminal coefficients.
     *
     *      Note that
     *          - i is in [0 ... N_P[       is the 'major' index (i.e. changes in i are not contiguous)
     *          - a is in [N_P ... K[       is the 'minor' index (i.e. changes in a are contiguous)
     */
    size_t vectorIndex(size_t i, size_t a) const;

    /**
     *  Calculate the Jacobian element with compound indices (i,a) and (k,c) at the given geminal coefficients @param g
     *
     *      i and k are subscripts, a and c are superscripts
     */
    double calculateJacobianElement(const Eigen::VectorXd& g, size_t i, size_t a, size_t k, size_t c) const;

    /**
     *  Calculate and return the Jacobian at the given geminal coefficients @param g
     */
    Eigen::MatrixXd calculateJacobian(const Eigen::VectorXd& g) const;

    /**
     *  Calculate the coordinate function at the given geminal coefficients @param g, with given indices.
     *
     *      i is the subscript and a is the superscript
     */
    double calculateCoordinateFunction(const Eigen::VectorXd& g, size_t i, size_t a) const;

    /**
     *  Calculate the coordinate functions for the pSEs at the given geminal coefficients @param g. This returns a vector F in which every entry is one of the coordinate functions
     */
    Eigen::VectorXd calculateCoordinateFunctions(const Eigen::VectorXd& g) const;

    /**
     *  Calculate the AP1roG energy given AP1roG geminal coefficients @param G
     */
    double calculateEnergy(const GQCG::AP1roGGeminalCoefficients& G) const;

    /**
     *  Set up and solve the projected Schrödinger equations for AP1roG
     */
    void solve();
};


}  // namespace GQCG



#endif /* AP1roGPSESolver_hpp */
