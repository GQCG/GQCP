#ifndef CartesianGTO_hpp
#define CartesianGTO_hpp


#include "math/ScalarFunction.hpp"
#include "math/LinearCombination.hpp"


namespace GQCP {

/**
 *  A class representing a Cartesian Gaussian-type orbital (GTO), which is often referred to as a 'primitive'
 *
 *  Mathematically speaking, a Cartesian GTO is a real-valued scalar function taking an Euclidean vector (3D-vector) as argument
 *
 *  Contracted GTOs can be expressed as linear combinations of GTOs: LinearCombination<CartesianGTO>
 */
class CartesianGTO : public ScalarFunction<double, double, 3> {
public:
    double alpha;  // exponent of the exponential
    double N;  // normalization factor
    std::array<size_t, 3> exponents;  // exponents of (x-X), (y-Y), (z-Z)
    Eigen::Vector3d center;  // center of the GTO (X, Y, Z)


public:
    // CONSTRUCTORS
    /**
     *  @param alpha        the exponent of the exponential
     *  @param exponents    the exponents of x, y and z, in that order
     *  @param center       the center of the Cartesian GTO
     */
    CartesianGTO(double alpha, const std::array<size_t, 3>& exponents, const Eigen::Vector3d& center);

    /**
     *  Default constructor setting everything to zero
     */
    CartesianGTO();


    // GETTERS
    double get_exponent() const { return this->alpha; }
    const std::array<size_t, 3>& get_exponents() const { return this->exponents; }
    const Eigen::Vector3d& get_center() const { return this->center; }


    // OPERATORS
    /**
     *  @param r        the value at which the GTO should be evaluated
     *
     *  @return the value of the GTO at the given position
     */
    double operator()(const Eigen::Vector3d& r) const override;


    // PUBLIC METHODS
    /**
     *  @param alpha   the exponent of the GTO
     *  @param c       the power of the Cartesian function x, y, z
     *
     *  @return one of the components of the total normalization factor
     */
    static double calculateNormalizationFactorComponent(double alpha, size_t c);

    /**
     *  @return the total normalization factor of the Cartesian GTO
     */
    double calculateNormalizationFactor() const;

    /**
     *  @param c        which component (x=0, y=1, z=2)
     *
     *  @return the derivative of this Cartesian GTO (with respect to the electronic coordinates) in the x-, y-, or z-direction
     */
    LinearCombination<CartesianGTO> calculateDerivative(size_t c) const;

    /**
     *  @return the gradient of this Cartesian GTO with respect to the electronic coordinates
     */
    Vector<LinearCombination<CartesianGTO>, 3> calculateGradient() const;
};


}  // namespace GQCP


#endif  /* CartesianGTO_hpp */
