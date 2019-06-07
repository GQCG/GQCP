#ifndef GQCP_NEWTONORBITALOPTIMIZER_HPP
#define GQCP_NEWTONORBITALOPTIMIZER_HPP


#include "OrbitalOptimization/BaseOrbitalOptimizer.hpp"
#include "OrbitalOptimization/OrbitalRotationGenerators.hpp"
#include "math/SquareMatrix.hpp"
#include "math/SquareRankFourTensor.hpp"

#include "Eigen/Dense"


namespace GQCP {


/**
 *  An intermediate abstract class that can be used to implement a Newton-step based orbital optimization.
 */
class NewtonOrbitalOptimizer : public BaseOrbitalOptimizer {
protected:
    VectorX<double> gradient;
    SquareMatrix<double> hessian;


public:
    // CONSTRUCTORS
    using BaseOrbitalOptimizer::BaseOrbitalOptimizer;  // inherit base constructors


    // PUBLIC PURE VIRTUAL METHODS

    /**
     *  @param ham_par      the current Hamiltonian parameters
     * 
     *  @return the current orbital gradient as a matrix
     */
    virtual SquareMatrix<double> calculateGradientMatrix(const HamiltonianParameters<double>& ham_par) const = 0;

    /**
     *  @param ham_par      the current Hamiltonian parameters
     * 
     *  @return the current orbital Hessian as a tensor
     */
    virtual SquareRankFourTensor<double> calculateHessianTensor(const HamiltonianParameters<double>& ham_par) const = 0;

    /**
     *  Use gradient and Hessian information to determine a new direction for the 'full' orbital rotation generators kappa. Note that a distinction is made between 'free' generators, i.e. those that are calculated from the gradient and Hessian information and the 'full' generators, which also include the redundant parameters (that can be set to zero). The 'full' generators are used to calculate the total rotation matrix using the matrix exponential
     * 
     *  @param ham_par      the current Hamiltonian parameters
     * 
     *  @return the new full set orbital generators, including the redundant parameters
     */
    virtual OrbitalRotationGenerators calculateNewFullOrbitalGenerators(const HamiltonianParameters<double>& ham_par) const = 0;


    // PUBLIC OVERRIDDEN METHODS

    /**
     *  Determine if the algorithm has converged or not
     *  Specifically for the Newton-step based algorithms, this function
     *      - computes the gradient and checks its norm for convergence
     *      - if the gradient is zero, the Hessian is calculated and diagonalized and positive/negative definiteness is checked
     * 
     *  @param ham_par      the current Hamiltonian parameters
     * 
     *  @return if the algorithm is considered to be converged
     */
    bool checkForConvergence(const HamiltonianParameters<double>& ham_par) override;

    /**
     *  Produce a new rotation matrix by either
     *      - continuing in the direction of the largest (in absolute value) non-conforming eigenvalue (i.e. the smallest (negative) eigenvalue for minimization algorithms and the largest (positive) eigenvalue for maximization algorithms)
     *      - using the Newton step if it is well-defined
     * 
     *  @param ham_par      the current Hamiltonian parameters
     * 
     *  @return a unitary matrix that will be used to rotate the current Hamiltonian parameters into the next iteration
     */
    SquareMatrix<double> calculateNewRotationMatrix(const HamiltonianParameters<double>& ham_par) override;


    // PUBLIC METHODS

    /**
     *  @param ham_par      the current Hamiltonian parameters
     * 
     *  @return the current orbital gradient as a vector. Matrix indices are converted to vector indices in the convention that p>q
     */
    VectorX<double> calculateGradientVector(const HamiltonianParameters<double>& ham_par) const;

    /**
     *  @param ham_par      the current Hamiltonian parameters
     * 
     *  @return the current orbital Hessian as a matrix
     */
    SquareMatrix<double> calculateHessianMatrix(const HamiltonianParameters<double>& ham_par) const;

    /**
     *  @return if a Newton step would be well-defined (i.e. the Hessian is positive definite for minimizations and negative definite for maximizations)
     */
    bool newtonStepIsWellDefined() const;

    /**
     *  If the Newton step is ill-defined, examine the Hessian and produce a new direction from it:
     *      - for minimization algorithms, this is the eigenvector that corresponds to the smallest (negative) eigenvalue of the Hessian
     *      - for maximization algorithms, this is the eigenvector that corresponds to the largest (positive) eigenvalue of the Hessian
     * 
     *  @return the new direction from the Hessian if the Newton step is ill-defined
     */
    VectorX<double> directionFromHessian() const;

    /**
     *  Use gradient and Hessian information to determine a new direction for the 'free' orbital rotation generators kappa. Note that a distinction is made between 'free' generators, i.e. those that are calculated from the gradient and Hessian information and the 'full' generators, which also include the redundant parameters (that can be set to zero). The 'full' generators are used to calculate the total rotation matrix using the matrix exponential
     * 
     *  @param ham_par      the current Hamiltonian parameters
     * 
     *  @return the new free orbital generators
     */
    OrbitalRotationGenerators calculateNewFreeOrbitalGenerators(const HamiltonianParameters<double>& ham_par) const;
};



}  // namespace GQCP



#endif  /* GQCP_NEWTONORBITALOPTIMIZER_HPP */
