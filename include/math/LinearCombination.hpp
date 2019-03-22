// This file is part of GQCG-gqcp.
// 
// Copyright (C) 2017-2019  the GQCG developers
// 
// GQCG-gqcp is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// GQCG-gqcp is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-gqcp.  If not, see <http://www.gnu.org/licenses/>.
// 
#ifndef LinearCombination_hpp
#define LinearCombination_hpp


#include "math/ScalarFunction.hpp"

#include <vector>


namespace GQCP {


/**
 *  A class template representing a linear combination of scalar functions
 *
 *  @tparam T                       the type of scalar function
 *  @tparam CoefficientScalar       the type of scalar that is used as coefficient
 */
template <typename CoefficientScalar, typename T>
class LinearCombination : public ScalarFunction<typename T::Valued, typename T::Scalar, T::Cols> {
    static_assert(std::is_base_of<ScalarFunction<typename T::Valued, typename T::Scalar, T::Cols>, T>::value, "LinearCombination: T must derive from ScalarFunction");


private:
    std::vector<CoefficientScalar> coefficients;
    std::vector<T> functions;


public:

    /*
     *  CONSTRUCTORS
     */

    /**
     *  @param coefficients     the coefficients of the linear combination
     *  @param functions        the scalar functions of the linear combination
     */
    LinearCombination(const std::vector<CoefficientScalar>& coefficients, const std::vector<T>& functions) :
        coefficients (coefficients),
        functions (functions)
    {
        if (coefficients.size() != functions.size()) {
            throw std::invalid_argument("LinearCombination(std::vector<CoefficientScalar>, std::vector<T>): the number of coefficients and functions should match");
        }
    }


    /**
     *  Default constructor: construct a 'zero vector', i.e. an empty linear combination
     */
    LinearCombination() :
        LinearCombination(std::vector<CoefficientScalar> {}, std::vector<T> {})
    {}


    /**
     *  Constructor providing a 'linear combination' of just one scalar function
     *
     *  @param coefficient      the one coefficient that belongs to the function
     *  @param function         one single scalar function
     */
    LinearCombination(CoefficientScalar coefficient, const T& function) :
        LinearCombination(std::vector<CoefficientScalar> {coefficient}, std::vector<T> {function})
    {}


    /**
     *  Constructor providing a 'linear combination' of just one scalar function, defaulting the coefficient to 1.0
     *
     *  @param function         one single scalar function
     */
    LinearCombination(const T& function) :
        LinearCombination(1.0, function)
    {}


    /**
     *  The constructor for a 'zero vector' given an integer argument
     *
     *  This constructor is added to fix errors in Eigen concerning 'Scalar(0)';
     */
    LinearCombination(int zero) :
        LinearCombination()
    {
        if (zero != 0) {
            throw std::invalid_argument("LinearCombination(int): Can't convert a non-zero integer to a zero vector");
        }
    }


    /*
     *  GETTERS
     */

    const std::vector<CoefficientScalar>& get_coefficients() const { return this->coefficients; }
    const std::vector<T>& get_functions() const { return this->functions; }


    /*
     *  OPERATORS implementing the notion of linear combinations
     */

    /**
     *  @param rhs      the right-hand side of the addition
     *
     *  @return the sum of this and the right-hand side
     */
    LinearCombination operator+(const LinearCombination& rhs) const {
        LinearCombination result = *this;
        result.append(rhs.coefficients, rhs.functions);
        return result;
    }

    /**
     *  @param rhs      the right-hand side of the addition
     *
     *  @return the sum of this and the right-hand side
     */
    LinearCombination& operator+=(const LinearCombination& rhs) {
        this->append(rhs.coefficients, rhs.functions);
        return *this;
    }

    /**
     *  @param scalar       the scalar to be used in the multiplication
     *
     *  @return the product of this with a scalar
     */
    LinearCombination operator*(CoefficientScalar scalar) const {

        if (std::abs(scalar) < 1.0e-16) {  // multiplication by zero returns a 'zero vector'
            return LinearCombination();
        }

        auto coefficients = this->coefficients;

        for (auto& coeff : coefficients) {
            coeff *= scalar;
        }

        return LinearCombination(coefficients, this->functions);
    }

    /**
     *  @param scalar       the scalar to be used in the multiplication
     *  @param rhs          the right-hand side of the multiplication
     *
     *  @return the product of a LinearCombination with a scalar
     */
    friend LinearCombination operator*(CoefficientScalar scalar, const LinearCombination& rhs) {
        return rhs.operator*(scalar);
    }

    using ScalarFunction<typename T::Valued, typename T::Scalar, T::Cols>::operator*;


    /**
     *  @return the negative of the linear combination
     */
    LinearCombination operator-() const {
        return this->operator*(-1.0);
    }

    /**
     *  @param rhs          the right-hand side of the subtraction
     *
     *  @return the difference of this and a right-hand side
     */
    LinearCombination operator-(const LinearCombination& rhs) const {
        return this->operator+(-rhs);
    }


    /*
     *  OTHER OPERATORS
     */

    /**
     *  @param x        the vector/point at which the scalar function is to be evaluated
     *
     *  @return the scalar function value of this linear combination at the given point
     */
    typename T::Valued operator()(const Vector<typename T::Scalar, T::Cols>& x) const override {
        size_t n = this->functions.size();

        CoefficientScalar value {};  // default initialization
        for (size_t i = 0; i < n; i++) {
            value += this->coefficients[i] * this->functions[i].operator()(x);  // evaluate every function of the linear combination
        }

        return value;
    }


    /**
     *  @param rhs      the right-hand side of the operator ==
     *
     *  @return whether two linear combinations are considered equal
     */
    bool operator==(const LinearCombination& rhs) const {
        return (this->coefficients == rhs.coefficients) && (this->functions == rhs.functions);
    }



    /*
     *  PUBLIC METHODS
     */
    
    /**
     *  Append the given coefficients and functions to this linear combination
     *
     *  @param coefficients     the coefficients that should be appended to this linear combination
     *  @param functions        the functions that should be appended to this linear combination
     */
    void append(const std::vector<CoefficientScalar>& coefficients, const std::vector<T>& functions) {

        if (coefficients.size() != functions.size()) {
            throw std::invalid_argument("LinearCombination::append(std::vector<CoefficientScalar>, std::vector<T>): the number of coefficients and functions should match");
        }

        this->coefficients.insert(this->coefficients.end(), coefficients.begin(), coefficients.end());
        this->functions.insert(this->functions.end(), functions.begin(), functions.end());
    }
};


}  // namespace GQCP



/*
 *  Make GQCP::LinearCombination<T> an Eigen scalar type
 */

namespace Eigen {

template<typename CoefficientScalar, typename T>
struct NumTraits<GQCP::LinearCombination<CoefficientScalar, T>> : public NumTraits<double> {  // permits to get the epsilon, dummy_precision, lowest, highest functions

    using Real = GQCP::LinearCombination<CoefficientScalar, T>;
    using NonInteger = GQCP::LinearCombination<CoefficientScalar, T>;
    using Nested = GQCP::LinearCombination<CoefficientScalar, T>;

    enum {
        IsComplex = 0,
        IsInteger = 0,
        IsSigned = 1,
        RequireInitialization = 1,
        ReadCost = 1,
        AddCost = 5,
        MulCost = 1000  // just put something big
    };
};


// Enable the scalar product of a LinearCombination<CoefficientScalar, T> with its own CoefficientScalar (both sides)

template<typename CoefficientScalar, typename T>
struct ScalarBinaryOpTraits<GQCP::LinearCombination<CoefficientScalar, T>, CoefficientScalar, Eigen::internal::scalar_product_op<GQCP::LinearCombination<CoefficientScalar, T>, CoefficientScalar>> {
    using ReturnType = GQCP::LinearCombination<CoefficientScalar, T>;
};

template<typename CoefficientScalar, typename T>
struct ScalarBinaryOpTraits<CoefficientScalar, GQCP::LinearCombination<CoefficientScalar, T>, Eigen::internal::scalar_product_op<CoefficientScalar, GQCP::LinearCombination<CoefficientScalar, T>>> {
    using ReturnType = GQCP::LinearCombination<CoefficientScalar, T>;
};


}  // namespace Eigen



#endif  /* LinearCombination_hpp */
