#ifndef GQCP_NUMERICALDERIVATIVE_HPP
#define GQCP_NUMERICALDERIVATIVE_HPP


#include "math.h"
#include "common.hpp"
#include "FockSpace/FockSpace.hpp"


namespace GQCP {


template<size_t T>
class NumericalDerivator {

private:
    double derivative;
    double function_value;
    NumericalDerivator<T-1> precursor;

public:
    NumericalDerivator (const UnaryFunction& uf, double start, double step_size) : precursor(uf, start, step_size){
        this->function_value = uf(start + step_size*T);
        this->derivative = 0;
        for (size_t i = 0; i < T+1; i++) {
            this->derivative += pow(-1,(i+1)) * FockSpace::calculateDimension(T,i) * this->get_function_value(i);
        }
        this->derivative /= pow(step_size, T);
        this->derivative *= pow(-1, T+1);
    }

    double get_derivative(size_t order = T) {
        if (order > T) {
            throw std::invalid_argument("requested derivative order was not computed");
        } else if  (order == T)  {
            return this->derivative;

        } else {
            return this->precursor.get_derivative(order);
        }
    }

    double get_function_value(size_t order = T) {
        if (order > T) {
            throw std::invalid_argument("function value of requested order was not computed");
        } else if  (order == T)  {
            return this->function_value;

        } else {
            return this->precursor.get_function_value(order);
        }
    }
};


template<>
class NumericalDerivator<0> {
private:
    double derivative;
    double function_value;

public:
    NumericalDerivator (const UnaryFunction& uf, double start, double step_size) {
        this->function_value = uf(start);
        this->derivative = this->function_value;
    }

    double get_derivative(size_t order = 0) {
        if (order != 0) {
            throw std::invalid_argument("requested derivative order was not computed");
        } else {
            return this->derivative;
        }
    }

    double get_function_value(size_t order = 0) {
        if (order != 0) {
            throw std::invalid_argument("function value of requested order was not computed");
        } else {
            return this->function_value;
        }
    }
};


}


#endif  // GQCP_NUMERICALDERIVATIVE_HPP
