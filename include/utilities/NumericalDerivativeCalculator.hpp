#ifndef GQCP_NUMERICALDERIVATIVE_HPP
#define GQCP_NUMERICALDERIVATIVE_HPP


#include "math.h"
#include "common.hpp"
#include "optimization/Eigenpair.hpp"
#include "FockSpace/FockSpace.hpp"
#include "miscellaneous.hpp"


namespace GQCP {

/**
 *  Class that recursively computes and stores T-order derivatives and the associated function values for unary functions
 */
template<size_t T>
class NumericalDerivativeCalculator {
protected:
    double derivative;
    double function_value;
    std::shared_ptr<NumericalDerivativeCalculator<T-1>> precursor;  // lower order DerivativeCalculator

    /**
     *  Numerically computes the derivative as : (-1)^(T+1) 1/(s^T) * sum^T_i=0 (-1)^(i+1) * (T)choose(i) * f(x + i*s)
     *
     *  @param s             a given step size to perform the T-order numerical derivative of a function
     *
     *  @return the T-th order derivative
     */
    double calculateDerivative(double s) {
        double derivative = 0;
        for (size_t i = 0; i < T+1; i++) {
            derivative += pow(-1,(i+1)) * binomialCoefficient(T,i) * this->get_function_value(i);
        }
        derivative /= pow(s, T);
        derivative *= pow(-1, T+1);
        return derivative;
    }

public:
    // CONSTRUCTORS
    /**
     *  Constructor for initialization of derived instances
     *
     * @param precursor
     */
    explicit NumericalDerivativeCalculator(std::shared_ptr<NumericalDerivativeCalculator<T-1>> precursor) : precursor(precursor) {};

    /**
     *  Recursively constructs T-th order down to 0-th order NumericalDerivativeCalculator
     *
     *  @param uf                    The Unary function we derive
     *  @param start                 starting parameter around which we derive
     *  @param step_size             step size for the numeric derivation approach
     */
    NumericalDerivativeCalculator (const UnaryFunction& uf, double start, double step_size) : precursor(std::make_shared<NumericalDerivativeCalculator<T-1>>(uf, start, step_size))
    {
        this->function_value = uf(start + step_size*T);
        this->derivative = calculateDerivative(step_size);
    }

    // GETTERS
    double get_derivative(size_t order = T) const {
        if (order > T) {
            throw std::invalid_argument("requested derivative order was not computed");
        } else if  (order == T)  {
            return this->derivative;

        } else {
            return this->precursor->get_derivative(order);
        }
    }

    double get_function_value(size_t order = T) const {
        if (order > T) {
            throw std::invalid_argument("function value of requested order was not computed");
        } else if  (order == T)  {
            return this->function_value;

        } else {
            return this->precursor->get_function_value(order);
        }
    }
};


/**
 *  Zero-th order template specialization
 */
template<>
class NumericalDerivativeCalculator<0> {
protected:
    double derivative;
    double function_value;

public:
    // CONSTRUCTORS
    NumericalDerivativeCalculator() = default;

    /**
     *  Constructs 0-th order NumericalDerivativeCalculator
     *
     *  @param uf                    The Unary function we derive
     *  @param start                 starting parameter around which we derive
     *  @param step_size             step size for the numeric derivation approach
     */
    NumericalDerivativeCalculator (const UnaryFunction& uf, double start, double step_size) {
        this->function_value = uf(start);
        this->derivative = this->function_value;
    }

    // GETTERS
    double get_derivative(size_t order = 0) const {
        if (order != 0) {
            throw std::invalid_argument("requested derivative order was not computed");
        } else {
            return this->derivative;
        }
    }

    double get_function_value(size_t order = 0) const {
        if (order != 0) {
            throw std::invalid_argument("function value of requested order was not computed");
        } else {
            return this->function_value;
        }
    }
};


/**
 *  Class that recursively computes and stores T-order derivatives and the associated function values and eigenvectors for numeric eigenproblems requiring a guess vector input
 */
using NumericEigenProblem = std::function<Eigenpair (double x, const Eigen::VectorXd& guess)>;

template<size_t T>
class NumericalGuessDerivativeCalculator : public NumericalDerivativeCalculator<T> {
private:
    Eigen::VectorXd eigenvector;

public:
    // CONSTRUCTORS
    /**
     *  Recursively constructs T-th order down to 0-th order NumericalGuessDerivativeCalculator
     *
     *  @param nef                   NumericEigenProblem function based on a single parameter
     *  @param start                 starting parameter around which we derive
     *  @param step_size             step size for the numeric derivation approach
     */
    NumericalGuessDerivativeCalculator (const NumericEigenProblem& uf, double start, double step_size, const Eigen::MatrixXd& guess) : NumericalDerivativeCalculator<T>(std::make_shared<NumericalGuessDerivativeCalculator<T-1>>(uf, start, step_size, guess))
    {
        const Eigenpair& eigenpair = uf(start + T*step_size, static_cast<const NumericalGuessDerivativeCalculator<T-1>&>(*this->precursor).get_eigenvector());
        this->function_value = eigenpair.get_eigenvalue();
        this->eigenvector = eigenpair.get_eigenvector();
        this->derivative = this->calculateDerivative(step_size);
    }

    // GETTERS
    const Eigen::VectorXd& get_eigenvector(size_t order = T) const {
        if (order > T) {
            throw std::invalid_argument("eigenvector of requested order was not computed");
        } else if  (order == T)  {
            return this->eigenvector;
        } else {
            return static_cast<const NumericalGuessDerivativeCalculator<T-1>&>(*this->precursor).get_eigenvector(order);
        }
    }
};


/**
 *  Zero-th order template specialization
 */
template<>
class NumericalGuessDerivativeCalculator<0> : public NumericalDerivativeCalculator<0> {
private:
    Eigen::VectorXd eigenvector;

public:
    // CONSTRUCTORS
    /**
     *  Constructs 0-th order NumericalGuessDerivativeCalculator
     *
     *  @param nef                   NumericEigenProblem function based on a single parameter
     *  @param start                 starting parameter around which we derive
     *  @param step_size             step size for the numeric derivation approach
     */
    NumericalGuessDerivativeCalculator (const NumericEigenProblem& nef, double start, double step_size, const Eigen::MatrixXd& guess) : NumericalDerivativeCalculator<0>() {
        const Eigenpair& eigenpair = nef(start, guess);
        this->function_value = eigenpair.get_eigenvalue();
        this->eigenvector = eigenpair.get_eigenvector();
        this->derivative = this->function_value;
    }

    // GETTERS
    const Eigen::VectorXd& get_eigenvector(size_t order = 0) const {
        if (order != 0) {
            throw std::invalid_argument("eigenvector of requested order was not computed");
        } else {
            return this->eigenvector;
        }
    }
};


}  // GQCP


#endif  // GQCP_NUMERICALDERIVATIVE_HPP
