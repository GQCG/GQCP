#include "AP1roG/AP1roGPSESolver.hpp"
#include <numopt.hpp>


namespace GQCG {


/*
 * CONSTRUCTORS
 */
/**
 *  Constructor based on a given number of electrons @param N and Hamiltonian parameters @param ham_par
 *
 *  The initial guess for the geminal coefficients is zero
 */
AP1roGPSESolver::AP1roGPSESolver(size_t N, const GQCG::HamiltonianParameters& ham_par) :
    K (ham_par.K),
    ham_par (ham_par),
    N_P (N / 2),
    initial_geminal_coefficients (GQCG::AP1roGGeminalCoefficients(this->N_P, this->K))
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
AP1roGPSESolver::AP1roGPSESolver(const GQCG::Molecule& molecule, const GQCG::HamiltonianParameters& ham_par) :
    AP1roGPSESolver(molecule.N, ham_par)
{}



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


/**
 *  Calculate the Jacobian element with compound indices (i,a) and (k,c) at the given geminal coefficients @param g
 *
 *      i and k are subscripts, a and c are superscripts
 */
double AP1roGPSESolver::calculateJacobianElement(const Eigen::VectorXd& g, size_t i, size_t a, size_t k, size_t c) const {

    Eigen::MatrixXd h_SO = this->ham_par.h.get_matrix_representation();
    Eigen::Tensor<double, 4> g_SO = this->ham_par.g.get_matrix_representation();

    double j_el = 0.0;


    // KISS implementation of the calculation of Jacobian elements
    if (i != k) {

        if (a != c) {  // i!=k and a!=c
            return 0.0;
        }

        else {  // i!=k and a == c
            j_el += g_SO(k,i,k,i) - g_SO(k,a,k,a) * g(this->vectorIndex(i, a));

            for (size_t b = this->N_P; b < this->K; b++) {
                if (b != a) {
                    j_el += g_SO(k,b,k,b) * g(this->vectorIndex(i, b));
                }
            }

        }
    }

    else {  // i==k

        if (a != c) {  // i==k and a!=c
            j_el += g_SO(a,c,a,c) - g_SO(i,c,i,c) * g(this->vectorIndex(i, a));

            for (size_t j = 0; j < this->N_P; j++) {
                if (j != i) {
                    j_el += g_SO(j,c,j,c) * g(this->vectorIndex(j, a));
                }
            }
        }

        else {  // i==k and a==c
            j_el += -2 * g_SO(a,i,a,i) * g(this->vectorIndex(i, a));

            for (size_t j = 0; j < this->N_P; j++) {
                if (j != i) {
                    j_el += 2 * (2 * g_SO(a,a,j,j) - g_SO(a,j,j,a)) - (2 * g_SO(i,i,j,j) - g_SO(i,j,j,i));
                }
            }

            j_el += 2 * (h_SO(a,a) - h_SO(i,i));

            j_el += g_SO(a,a,a,a) - g_SO(i,i,i,i);

            for (size_t b = this->N_P; b < this->K; b++) {
                if (b != a) {
                    j_el += - g_SO(i,b,i,b) * g(this->vectorIndex(i, b));
                }
            }

            for (size_t j = 0; j < this->N_P; j++) {
                if (j != i) {
                    j_el += - g_SO(j,a,j,a) * g(this->vectorIndex(j, a));
                }
            }
        }

    }

    return j_el;
}


/**
 *  Calculate and return the Jacobian at the given geminal coefficients @param g
 */
Eigen::MatrixXd AP1roGPSESolver::calculateJacobian(const Eigen::VectorXd& g) const {

    // The dimension of the Jacobian is (K-N_P)*N_P times (K-N_P)*N_P
    Eigen::MatrixXd J = Eigen::MatrixXd::Zero((this->K - this->N_P) * this->N_P, (this->K - this->N_P) * this->N_P);

    // Loop over all Jacobian elements to construct it
    size_t number_of_geminal_coefficients = AP1roGGeminalCoefficients::numberOfGeminalCoefficients(N_P, K);
    for (size_t mu = 0; mu < number_of_geminal_coefficients; mu++) {
        for (size_t nu = 0; nu < number_of_geminal_coefficients; nu++) {

            // Convert the vector indices mu and nu into matrix indices
            size_t i = this->matrixIndexMajor(nu);
            size_t a = this->matrixIndexMinor(nu);
            size_t k = this->matrixIndexMajor(mu);
            size_t c = this->matrixIndexMinor(mu);

            J(mu, nu) = this->calculateJacobianElement(g, i, a, k, c);
        }
    }

    return J;
}


/**
 *  Calculate the coordinate function at the given geminal coefficients @param g, with given indices.
 *
 *      i is the subscript and a is the superscript
 */
double AP1roGPSESolver::calculateCoordinateFunction(const Eigen::VectorXd& g, size_t i, size_t a) const {

    Eigen::MatrixXd h_SO = this->ham_par.h.get_matrix_representation();
    Eigen::Tensor<double, 4> g_SO = this->ham_par.g.get_matrix_representation();

    double f = 0.0;

    // A KISS implementation of the AP1roG pSE equations
    f += g_SO(a,i,a,i) * (1 - std::pow(g(this->vectorIndex(i, a)), 2));

    for (size_t j = 0; j < this->N_P; j++) {
        if (j != i) {
            f += 2 * ((2 * g_SO(a,a,j,j) - g_SO(a,j,j,a)) - (2 * g_SO(i,i,j,j) - g_SO(i,j,j,i))) * g(this->vectorIndex(i, a));
        }
    }

    f += 2 * (h_SO(a,a) - h_SO(i,i)) * g(this->vectorIndex(i, a));

    f += (g_SO(a,a,a,a) - g_SO(i,i,i,i)) * g(this->vectorIndex(i, a));

    for (size_t b = this->N_P; b < this->K; b++) {
        if (b != a) {
            f += (g_SO(a,b,a,b) - g_SO(i,b,i,b) * g(this->vectorIndex(i, a))) * g(this->vectorIndex(i, b));
        }
    }

    for (size_t j = 0; j < this->N_P; j++) {
        if (j != i) {
            f += (g_SO(j,i,j,i) - g_SO(j,a,j,a) * g(this->vectorIndex(i, a))) * g(this->vectorIndex(j, a));
        }
    }

    for (size_t b = this->N_P; b < this->K; b++) {
        if (b != a) {

            for (size_t j = 0; j < this->N_P; j++) {
                if (j != i) {
                    f += g_SO(j,b,j,b) * g(this->vectorIndex(j, a)) * g(this->vectorIndex(i, b));
                }
            }

        }
    }

    return f;
}


/**
 *  Calculate the coordinate functions for the pSEs at the given geminal coefficients @param g. This returns a vector F in which every entry is one of the coordinate functions
 */
Eigen::VectorXd AP1roGPSESolver::calculateCoordinateFunctions(const Eigen::VectorXd& g) const {

    GQCG::AP1roGGeminalCoefficients gem_coeff (g, this->N_P, this->K);
    return this->calculateCoordinateFunctions(gem_coeff);
}


/**
 *  Calculate the coordinate functions or the PSEs at the given geminal coefficients @param G. @returns a vector F in which every entry is one of the coordinate functions
 */
Eigen::VectorXd AP1roGPSESolver::calculateCoordinateFunctions(const GQCG::AP1roGGeminalCoefficients& G) const {

    // The dimension of the vector F is the number of coordinate functions, which is (K-N_P)*N_P
    Eigen::VectorXd F = Eigen::VectorXd::Zero ((this->K - this->N_P) * this->N_P);

    Eigen::VectorXd g = G.asVector();

    // Loop over all the F elements to construct it
    size_t number_of_geminal_coefficients = AP1roGGeminalCoefficients::numberOfGeminalCoefficients(N_P, K);
    for (size_t mu = 0; mu < number_of_geminal_coefficients; mu++) {

        // Convert the vector indices mu into matrix indices
        size_t i = this->matrixIndexMajor(mu);
        size_t a = this->matrixIndexMinor(mu);

        F(mu) = this->calculateCoordinateFunction(g, i, a);
    }

    return F;
}



/**
 *  Calculate the AP1roG energy given AP1roG geminal coefficients @param G
 */
double AP1roGPSESolver::calculateEnergy(const GQCG::AP1roGGeminalCoefficients& G) const {

    Eigen::MatrixXd h_SO = this->ham_par.h.get_matrix_representation();
    Eigen::Tensor<double, 4> g_SO = this->ham_par.g.get_matrix_representation();
    Eigen::VectorXd g = G.asVector();

    double E = 0.0;

    // KISS implementation of the AP1roG energy
    for (size_t j = 0; j < this->N_P; j++) {
        E += 2 * h_SO(j,j);

        for (size_t k = 0; k < this->N_P; k++) {
            E += 2 * g_SO(k,k,j,j) - g_SO(k,j,j,k);
        }

        for (size_t b = this->N_P; b < this->K; b++) {
            E += g_SO(j,b,j,b) * g(this->vectorIndex(j, b));
        }
    }

    return E;
}


/**
 *  Set up and solve the projected Schrödinger equations for AP1roG
 */
void AP1roGPSESolver::solve() {

    // Solve the AP1roG equations using a Newton-based algorithm


    // using std::function: see (http://en.cppreference.com/w/cpp/utility/functional/function)
    // using Lambdas: see (https://stackoverflow.com/a/44792017)
    numopt::VectorFunction f = [this](const Eigen::VectorXd& x) {return this->calculateCoordinateFunctions(x); };
    numopt::JacobianFunction J = [this](const Eigen::VectorXd& x) {return this->calculateJacobian(x); };


    Eigen::VectorXd x0 = this->initial_geminal_coefficients.asVector();
    numopt::syseq::NewtonSystemOfEquationsSolver syseq_solver (x0, f, J);
    syseq_solver.solve();


    // Set the solution
    GQCG::AP1roGGeminalCoefficients geminal_coefficients (syseq_solver.get_solution(), this->N_P, this->K);
    double electronic_energy = this->calculateEnergy(geminal_coefficients);
    this->solution = GQCG::AP1roG(geminal_coefficients, electronic_energy);
}


}  // namespace GQCG
