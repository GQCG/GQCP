#include "AP1roG/AP1roGJacobiOrbitalOptimizer.hpp"

namespace GQCG {


/*
 *  CON
 */
/*
 * CONSTRUCTORS
 */
/**
 *  Constructor based on a given number of electrons @param N and Hamiltonian parameters @param ham_par
 *
 *  The initial guess for the geminal coefficients is zero
 */
AP1roGJacobiOrbitalOptimizer::AP1roGJacobiOrbitalOptimizer(size_t N, const GQCG::HamiltonianParameters& ham_par, double oo_threshold) :
    K (ham_par.K),
    ham_par (ham_par),
    N_P (N / 2),
    oo_threshold (oo_threshold)
{
    // Check if we have an even number of electrons
    if ((N % 2) != 0) {
        throw std::invalid_argument("The given number of electrons is odd.");
    }
}


/**
 *  Constructor based on a given @param molecule and Hamiltonian parameters @param ham_par
 *
 *  The initial guess for the geminal coefficients is zero
 */
AP1roGJacobiOrbitalOptimizer::AP1roGJacobiOrbitalOptimizer(const GQCG::Molecule& molecule, const GQCG::HamiltonianParameters& ham_par, double oo_threshold) :
    AP1roGJacobiOrbitalOptimizer(molecule.N, ham_par, oo_threshold)
{}



/*
 *  PUBLIC METHODS
 */
/**
 *  Given the two indices of spatial orbitals @param p and @param q that will be Jacobi-rotated, calculate the
 *  coefficients (which are @members)
 *      - A1, B1, C1            to be used in occupied-occupied rotations
 *      - A2, B2, C2, D2, E2    to be used in occupied-virtual rotations
 *      - A3, B3, C3            to be used in virtual-virtual rotations
 */
//void AP1roGJacobiOrbitalOptimizer::calculateJacobiCoefficients(size_t p, size_t q) {
//
//    // Implementation of the Jacobi rotation coefficients with disjoint cases for p and q
//
//    // Occupied-occupied rotations: if p <= N_P and q <= N_P for computers
//    if ((p < this->N_P) && (q < this->N_P)) {
//
//        this->A1 = 0.0;
//        this->B1 = 0.0;
//        this->C1 = 0.0;
//
//
//        for (size_t b = this->N_P; b < this->K; b++) {
//            this->A1 -= 0.5 * (this->so_basis.get_g_SO(b,p,b,p) - this->so_basis.get_g_SO(b,q,b,q)) * (this->g(this->vectorIndex(p, b)) - this->g(this->vectorIndex(q, b)));
//            this->B1 += 0.5 * (this->so_basis.get_g_SO(b,p,b,p) - this->so_basis.get_g_SO(b,q,b,q)) * (this->g(this->vectorIndex(p, b)) - this->g(this->vectorIndex(q, b)));
//            this->C1 += this->so_basis.get_g_SO(b,p,b,q) * (this->g(this->vectorIndex(q, b)) - this->g(this->vectorIndex(p, b)));
//        }
//    }
//
//
//    // Occupied-virtual rotations: if p > N_P and q <= N_P for computers
//    else if ((p >= this->N_P) && (q < this->N_P)) {
//
//        this->A2 = 0.0;
//        this->B2 = 0.0;
//        this->C2 = 0.0;
//        this->D2 = 0.0;
//        this->E2 = 0.0;
//
//
//        this->A2 += this->so_basis.get_h_SO(p,p) - this->so_basis.get_h_SO(q,q) + 0.375 * (this->so_basis.get_g_SO(p,p,p,p) + this->so_basis.get_g_SO(q,q,q,q)) * (1 - this->g(this->vectorIndex(q, p))) - 0.25 * this->so_basis.get_g_SO(p,p,q,q) * (7 + this->g(this->vectorIndex(q, p))) + 0.5 * this->so_basis.get_g_SO(p,q,p,q) * (3 + this->g(this->vectorIndex(q, p)));
//        this->B2 += this->so_basis.get_h_SO(q,q) - this->so_basis.get_h_SO(p,p) + 2 * this->so_basis.get_g_SO(p,p,q,q) + 0.5 * (this->so_basis.get_g_SO(p,p,p,p) + this->so_basis.get_g_SO(q,q,q,q)) * (this->g(this->vectorIndex(q, p)) - 1) - this->so_basis.get_g_SO(p,q,p,q) * (1 + this->g(this->vectorIndex(q, p)));
//        this->C2 += 2 * this->so_basis.get_h_SO(p,q) + (this->so_basis.get_g_SO(p,p,p,q) - this->so_basis.get_g_SO(p,q,q,q)) * (1 - this->g(this->vectorIndex(q, p)));
//        this->D2 += 0.125 * (this->so_basis.get_g_SO(p,p,p,p) + this->so_basis.get_g_SO(q,q,q,q) - 2 * (this->so_basis.get_g_SO(p,p,q,q) + 2 * this->so_basis.get_g_SO(p,q,p,q))) * (1 - this->g(this->vectorIndex(q, p)));
//        this->E2 += 0.5 * (this->so_basis.get_g_SO(p,p,p,q) - this->so_basis.get_g_SO(p,q,q,q)) * (this->g(this->vectorIndex(q, p)) - 1);
//
//        for (size_t j = 0; j < this->N_P; j++) {
//            this->A2 += 2 * (this->so_basis.get_g_SO(j,j,p,p) - this->so_basis.get_g_SO(j,j,q,q)) - 0.5 * (this->so_basis.get_g_SO(j,p,j,p) - this->so_basis.get_g_SO(j,q,j,q)) * (2 + this->g(this->vectorIndex(j, p)));
//            this->B2 += 2 * (this->so_basis.get_g_SO(j,j,q,q) - this->so_basis.get_g_SO(j,j,p,p)) + 0.5 * (this->so_basis.get_g_SO(j,p,j,p) - this->so_basis.get_g_SO(j,q,j,q)) * (2 + this->g(this->vectorIndex(j, p)));
//            this->C2 += 4 * this->so_basis.get_g_SO(j,j,p,q) - this->so_basis.get_g_SO(j,p,j,q) * (2 + this->g(this->vectorIndex(j, p)));
//        }
//
//        for (size_t b = this->N_P; b < this->K; b++) {
//            this->A2 += 0.5 * (this->so_basis.get_g_SO(b,p,b,p) - this->so_basis.get_g_SO(b,q,b,q)) * this->g(this->vectorIndex(q, b));
//            this->B2 += 0.5 * (this->so_basis.get_g_SO(b,q,b,q) - this->so_basis.get_g_SO(b,p,b,p)) * this->g(this->vectorIndex(q, b));
//            this->C2 += this->so_basis.get_g_SO(b,p,b,q) * this->g(this->vectorIndex(q, b));
//        }
//    }
//
//
//    // Virtual-virtual rotations: if p > N_P and q > N_P for computers
//    else if ((p >= this->N_P) && (q >= this->N_P )) {
//
//        this->A3 = 0.0;
//        this->B3 = 0.0;
//        this->C3 = 0.0;
//
//
//        for (size_t j = 0; j < this->N_P; j++) {
//            this->A3 -= 0.5 * (this->so_basis.get_g_SO(j,p,j,p) - this->so_basis.get_g_SO(j,q,j,q)) * (this->g(this->vectorIndex(j, p)) - this->g(this->vectorIndex(j, q)));
//            this->B3 += 0.5 * (this->so_basis.get_g_SO(j,p,j,p) - this->so_basis.get_g_SO(j,q,j,q)) * (this->g(this->vectorIndex(j, p)) - this->g(this->vectorIndex(j, q)));
//            this->C3 += this->so_basis.get_g_SO(j,p,j,q) * (this->g(this->vectorIndex(j, q)) - this->g(this->vectorIndex(j, p)));
//        }
//    }
//
//
//    else {  // this means that p <= q
//        throw std::invalid_argument("The given p and q are invalid: p must be larger than q.");
//    }
//
//    this->are_calculated_jacobi_coefficients = true;
//}


}  // namespace GQCG
