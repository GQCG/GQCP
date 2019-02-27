#ifndef LinearCombination_hpp
#define LinearCombination_hpp


#include "math/MultipliableScalarFunction.hpp"



namespace GQCP {


/**
 *  @tparam T    a type of scalar function
 *
 *  TODO: T should be a scalar function?
 *  TODO: double should also become a template?
 */
template <typename T>
class LinearCombination : public MultipliableScalarFunction<typename T::Valued, typename T::Scalar, T::Cols> {
public:
    std::vector<double> coefficients;
    std::vector<T> functions;


    // CONSTRUCTORS
    LinearCombination(const std::vector<double>& coefficients, const std::vector<T>& functions) :
        coefficients (coefficients),
        functions (functions)
    {}


    /**
     *  Default constructor: construct a zero vector
     */
    LinearCombination() :
        LinearCombination(std::vector<double> {0.0}, std::vector<T> {T()})
    {}


    /**
     *
     */
    LinearCombination(double coefficient, const T& function) :
        LinearCombination(std::vector<double> {coefficient}, std::vector<T> {function})
    {}


    /**
     *  Constructor providing a 'linear combination' of just one scalar function
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
            throw std::invalid_argument("Can't convert a non-zero integer to a zero vector");
        }
    }


    const std::vector<double>& get_coefficients() const { return this->coefficients; }
    const std::vector<T>& get_functions() const { return this->functions; }


    /*
     *  OPERATORS implementing the notion of linear combinations
     */

    /**
     *  Addition
     */
    LinearCombination operator+(const LinearCombination& rhs) const {
        LinearCombination result = *this;
        result.append(rhs.coefficients, rhs.functions);
        return result;
    }

    LinearCombination& operator+=(const LinearCombination& rhs) {
        this->append(rhs.coefficients, rhs.functions);
        return *this;
    }

    /**
     *  Scalar multiplication
     */
    LinearCombination operator*(double scalar) const {
        auto coefficients = this->coefficients;

        for (auto& coeff : coefficients) {
            coeff *= scalar;
        }

        return LinearCombination(coefficients, this->functions);
    }

    friend LinearCombination operator*(double scalar, const LinearCombination& rhs) {
        return rhs.operator*(scalar);
    }

    using MultipliableScalarFunction<typename T::Valued, typename T::Scalar, T::Cols>::operator*;


    /**
     *  Negation as scalar multiplication by -1
     */
    LinearCombination operator-() const {
        return this->operator*(-1);
    }

    /**
     *  Subtraction as negated addition
     */
    LinearCombination operator-(const LinearCombination& rhs) const {
        return this->operator+(-rhs);
    }



    // OPERATOR() implements the notion of a scalar function
    typename T::Valued operator()(const Eigen::Matrix<typename T::Scalar, T::Cols, 1>& x) const override {
        size_t n = this->functions.size();

        double value = 0.0;
        for (size_t i = 0; i < n; i++) {
            value += this->coefficients[i] * this->functions[i].operator()(x);  // evaluate every function of the linear combination
        }

        return value;
    }



    /**
     *  Overloading of operator<< to work with LinearCombinations
     */
    friend std::ostream& operator<<(std::ostream& os, const LinearCombination<T>& lc) {
        size_t n = lc.functions.size();

        for (size_t i = 0; i < n; i++) {
            os << lc.coefficients[i] << "*(" << lc.functions[i] << ')';

            if (i != n-1) {
                os << " + ";
            }
        }

        return os;
    }


    // PUBLIC METHODS
    /**
     *  Append the given coefficients and functions to this linear combination
     */
    void append(const std::vector<double>& coefficients, const std::vector<T>& functions) {
        this->coefficients.insert(this->coefficients.end(), coefficients.begin(), coefficients.end());
        this->functions.insert(this->functions.end(), functions.begin(), functions.end());
    }
};


}  // namespace GQCP




/*
 *  Make GQCP::LinearCombination<T> an Eigen scalar type
 */

namespace Eigen {

template<typename T>
struct NumTraits<GQCP::LinearCombination<T>> : public NumTraits<double> {  // permits to get the epsilon, dummy_precision, lowest, highest functions

    using Real = GQCP::LinearCombination<T>;
    using NonInteger = GQCP::LinearCombination<T>;
    using Nested = GQCP::LinearCombination<T>;

    enum {
        IsComplex = 0,
        IsInteger = 0,
        IsSigned = 1,
        RequireInitialization = 1,
        ReadCost = 1,
        AddCost = 5,
        MulCost = 1000
    };
};


// Enable custom scalar * double
template<typename T>
struct Eigen::ScalarBinaryOpTraits<GQCP::LinearCombination<T>, double, Eigen::internal::scalar_product_op<GQCP::LinearCombination<T>, double>> {
    using ReturnType = GQCP::LinearCombination<T>;
};

// Enable double * custom scalar
template<typename T>
struct Eigen::ScalarBinaryOpTraits<double, GQCP::LinearCombination<T>, Eigen::internal::scalar_product_op<double, GQCP::LinearCombination<T>>> {
    using ReturnType = GQCP::LinearCombination<T>;
};


}  // namespace Eigen



#endif  /* LinearCombination_hpp */
