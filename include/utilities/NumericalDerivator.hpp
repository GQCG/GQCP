#ifndef GQCP_NUMERICALDERIVATIVE_HPP
#define GQCP_NUMERICALDERIVATIVE_HPP


#include "math.h"
#include "common.hpp"
#include "optimization/Eigenpair.hpp"
#include "FockSpace/FockSpace.hpp"


namespace GQCP {

/**
 *  Class that recursively computes and stores T-order derivatives and the associated function values for unary functions
 */
template<size_t T>
class NumericalDerivator {
protected:
    double derivative;
    double function_value;
    std::shared_ptr<NumericalDerivator<T-1>> precursor;

    /**
     *  Numerically computes the derivative as : (-1)^(T+1) 1/(s^T) * sum^T_i=0 (-1)^(i+1) * (T)choose(i) * f(x + i*s)
     *
     *  @param s             a given step size to perform the T-order numerical derivative of a function
     *  @return the T-th order derivative
     */
    double calculateDerivative(double s) {
        double derivative = 0;
        for (size_t i = 0; i < T+1; i++) {
            derivative += pow(-1,(i+1)) * FockSpace::calculateDimension(T,i) * this->get_function_value(i);
        }
        derivative /= pow(s, T);
        derivative *= pow(-1, T+1);
        return derivative;
    }

public:
    explicit NumericalDerivator(std::shared_ptr<NumericalDerivator<T-1>> precursor) : precursor(precursor) {};

    NumericalDerivator (const UnaryFunction& uf, double start, double step_size) : precursor(std::make_shared<NumericalDerivator<T-1>>(uf, start, step_size)){
        this->function_value = uf(start + step_size*T);
        this->derivative = calculateDerivative(step_size);
    }

    /**
     *  @param order                order of the derivative that is requested
     *  @return the derivative of the requested order
     */
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
class NumericalDerivator<0> {
protected:
    double derivative;
    double function_value;

public:
    NumericalDerivator() = default;

    NumericalDerivator (const UnaryFunction& uf, double start, double step_size) {
        this->function_value = uf(start);
        this->derivative = this->function_value;
    }

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
 *  Class that recursively computes and stores T-order derivatives and the associated function values and eigenvectors for unary numeric eigen problems requiring a guess vector input
 */
using NumericEigenProblem = std::function<Eigenpair (double x, const Eigen::VectorXd& guess)>;

template<size_t T>
class NumericalGuessDerivator : public NumericalDerivator<T> {
private:
    Eigen::VectorXd eigenvector;

public:
    NumericalGuessDerivator (const NumericEigenProblem& uf, double start, double step_size, const Eigen::MatrixXd& guess) :
    NumericalDerivator<T>(std::make_shared<NumericalGuessDerivator<T-1>>(uf, start, step_size, guess))
    {
        const Eigenpair& eigenpair = uf(start + T*step_size, static_cast<const NumericalGuessDerivator<T-1>&>(*this->precursor).get_eigenvector());
        this->function_value = eigenpair.get_eigenvalue();
        this->eigenvector = eigenpair.get_eigenvector();
        this->derivative = this->calculateDerivative(step_size);
    }

    const Eigen::VectorXd& get_eigenvector(size_t order = T) const {
        if (order > T) {
            throw std::invalid_argument("eigenvector of requested order was not computed");
        } else if  (order == T)  {
            return this->eigenvector;
        } else {
            return static_cast<const NumericalGuessDerivator<T-1>&>(*this->precursor).get_eigenvector(order);
        }
    }
};

/**
 *  Zero-th order template specialization
 */
template<>
class NumericalGuessDerivator<0> : public NumericalDerivator<0> {
    Eigen::VectorXd eigenvector;
public:
    NumericalGuessDerivator (const NumericEigenProblem& uf, double start, double step_size, const Eigen::MatrixXd& guess) : NumericalDerivator<0>() {
        const Eigenpair& eigenpair = uf(start, guess);
        this->function_value = eigenpair.get_eigenvalue();
        this->eigenvector = eigenpair.get_eigenvector();
        this->derivative = this->function_value;
    }

    const Eigen::VectorXd& get_eigenvector(size_t order = 0) const {
        if (order != 0) {
            throw std::invalid_argument("eigenvector of requested order was not computed");
        } else {
            return this->eigenvector;
        }
    }

};






}


#endif  // GQCP_NUMERICALDERIVATIVE_HPP
