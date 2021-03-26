// This file is part of GQCG-GQCP.
//
// Copyright (C) 2017-2020  the GQCG developers
//
// GQCG-GQCP is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// GQCG-GQCP is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-GQCP.  If not, see <http://www.gnu.org/licenses/>.

#pragma once


#include "Mathematical/Functions/Function.hpp"
#include "Mathematical/Functions/VectorSpaceArithmetic.hpp"

#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>

#include <vector>


namespace GQCP {


/**
 *  A type that represents a linear combination of (homogeneous) mathematical functions.
 *
 *  @tparam _Coefficient        The type of a coefficient.
 *  @tparam _FunctionType       The type of the function.
 * 
 *  @note The type of the functions must derive from GQCP::Function, i.e. it must be evaluable.
 */
template <typename _Coefficient, typename _FunctionType>
class EvaluableLinearCombination:
    public Function<typename _FunctionType::OutputType, typename _FunctionType::InputType>,
    public VectorSpaceArithmetic<EvaluableLinearCombination<_Coefficient, _FunctionType>, _Coefficient> {

public:
    // The type of a coefficient.
    using Coefficient = _Coefficient;

    // The type of the function.
    using FunctionType = _FunctionType;

    // The return type of the `operator()`.
    using OutputType = typename FunctionType::OutputType;

    // The input type of the `operator()`.
    using InputType = typename FunctionType::InputType;

    static_assert(std::is_base_of<Function<OutputType, InputType>, FunctionType>::value, "EvaluableLinearCombination: FunctionType must derive from `Function`.");

    // The type of 'this'.
    using Self = EvaluableLinearCombination<Coefficient, FunctionType>;

protected:
    // The coefficients of the linear combination.
    std::vector<Coefficient> m_coefficients;

    // The functions of the linear combination.
    std::vector<FunctionType> m_functions;


public:
    /*
     *  MARK: Constructors
     */

    /**
     *  @param coefficients         The coefficients of the linear combination.
     *  @param functions            The functions of the linear combination.
     */
    EvaluableLinearCombination(const std::vector<Coefficient>& coefficients, const std::vector<FunctionType>& functions) :
        m_coefficients {coefficients},
        m_functions {functions} {

        if (coefficients.size() != functions.size()) {
            throw std::invalid_argument("EvaluableLinearCombination(const std::vector<Coefficient>&, const std::vector<FunctionType>&): The number of coefficients and functions do not match.");
        }
    }


    /**
     *  The default constructor. Construct a 'zero vector', i.e. an empty linear combination.
     */
    EvaluableLinearCombination() :
        EvaluableLinearCombination(std::vector<Coefficient> {}, std::vector<FunctionType> {}) {}


    /**
     *  Construct a 'linear combination' of just one function.
     *
     *  @param coefficient          The one coefficient that belongs to the function.
     *  @param function             The one single function.
     */
    EvaluableLinearCombination(const Coefficient coefficient, const FunctionType& function) :
        EvaluableLinearCombination(std::vector<Coefficient> {coefficient}, std::vector<FunctionType> {function}) {}


    /**
     *  Construct a 'linear combination' of just one function, setting the coefficient to 1.0.
     *
     *  @param function             The one single function.
     */
    EvaluableLinearCombination(const FunctionType& function) :
        EvaluableLinearCombination(Coefficient {1.0}, function) {}


    /**
     *  Construct a 'zero vector' given the '0' integer literal.
     * 
     *  @param zero         The '0' integer literal.
     */
    EvaluableLinearCombination(const int zero) :
        EvaluableLinearCombination() {

        if (zero != 0) {
            throw std::invalid_argument("EvaluableLinearCombination(const int): Can't convert a non-zero integer to a 'zero' instance.");
        }
    }


    /*
     *  MARK: Vector space arithmetic
     */

    /**
     *  Addition-assignment.
     */
    Self& operator+=(const Self& rhs) override {
        this->append(rhs.m_coefficients, rhs.m_functions);
        return *this;
    }


    /**
     *  Scalar multiplication-assignment.
     */
    Self& operator*=(const Coefficient& a) override {

        if (std::abs(a) < 1.0e-16) {  // multiplication by zero returns a 'zero vector'
            this->m_coefficients = std::vector<Coefficient> {};
            this->m_functions = std::vector<FunctionType> {};
            return *this;
        }

        std::transform(this->m_coefficients.begin(), this->m_coefficients.end(),
                       this->m_coefficients.begin(), [a](const Coefficient& C) { return C * a; });

        return *this;
    }


    /*
     *  MARK: General information
     */

    /**
     *  @return A textual description of this linear combination.
     */
    std::string description() const {

        std::string description = "[";
        for (size_t i = 0; i < this->length(); i++) {

            // Provide the coefficient and the function's description.
            description += (boost::format("(%|.3f|) %s") % this->coefficient(i) % this->function(i).description()).str();

            // Add a '+' when we're not at the end.
            if (i != this->length() - 1) {
                description += " + ";
            } else {
                description += "]";
            }
        }

        return description;
    }

    /**
     *  @return The length of the linear combination, i.e. the number of coefficients/functions inside it.
     */
    size_t length() const { return this->m_coefficients.size(); }


    /*
     *  MARK: Access
     */

    /**
     *  Access a coefficient.
     * 
     *  @param i        An index.
     * 
     *  @return A read-only reference to the coefficient of this linear combination that corresponds to the given index.
     */
    const Coefficient& coefficient(const size_t i) const { return this->m_coefficients[i]; }

    /**
     *  @return All the coefficients of the linear combination.
     */
    const std::vector<Coefficient>& coefficients() const { return this->m_coefficients; }

    /**
     *  Access a function.
     * 
     *  @param i        An index.
     * 
     *  @return A read-only reference to the function of this linear combination that corresponds to the given index.
     */
    const FunctionType& function(const size_t i) const { return this->m_functions[i]; }

    /**
     *  @return All the scalar functions of this linear combination.
     */
    const std::vector<FunctionType>& functions() const { return this->m_functions; }


    /*
     *  MARK: Conforming to `Function`
     */

    /**
     *  Evaluate this linear combination.
     * 
     *  @param in           The argument for which the function is to be evaluated.
     *
     *  @return The function value of this linear combination for the given argument.
     */
    OutputType operator()(const InputType& in) const override {

        // Evaluate every function of the linear combination.
        OutputType out {};  // Default initialization of the `OutputType`.
        const auto n = this->length();
        for (size_t i = 0; i < n; i++) {
            out += this->m_coefficients[i] * this->m_functions[i].operator()(in);
        }

        return out;
    }


    /*
     *  MARK: Comparison operators
     */

    /**
     *  @param other        The linear combination that this one is compared to.
     *
     *  @return Whether this linear combination is equal to the other.
     */
    bool operator==(const Self& other) const {
        return (this->m_coefficients == other.m_coefficients) && (this->m_functions == other.m_functions);
    }


    /*
     *  MARK: Appending
     */

    /**
     *  Append the given coefficient and function to this linear combination.
     *
     *  @param coefficients         The coefficient that should be appended to this linear combination.
     *  @param functions            The function that should be appended to this linear combination.
     *  @param threshold            The threshold for the (absolute value of the) coefficient in order to be included in this linear combination.
     */
    void append(const Coefficient& coefficient, const FunctionType& function, const double threshold = 1.0e-16) {

        // Only enlarge the linear combination for sufficiently large coefficients.
        if (std::abs(coefficient) < threshold) {
            return;
        } else {
            this->m_coefficients.push_back(coefficient);
            this->m_functions.push_back(function);
        }
    }


    /**
     *  Append the given coefficients and functions to this linear combination.
     *
     *  @param coefficients         The coefficients that should be appended to this linear combination.
     *  @param functions            The functions that should be appended to this linear combination.
     */
    void append(const std::vector<Coefficient>& coefficients, const std::vector<FunctionType>& functions) {

        if (coefficients.size() != functions.size()) {
            throw std::invalid_argument("EvaluableLinearCombination::append(const std::vector<Coefficient>&, const std::vector<FunctionType>&): The number of coefficients and functions do not match.");
        }

        this->m_coefficients.insert(this->m_coefficients.end(), coefficients.begin(), coefficients.end());
        this->m_functions.insert(this->m_functions.end(), functions.begin(), functions.end());
    }
};


}  // namespace GQCP


#include <Eigen/Dense>

/*
 *  Make GQCP::EvaluableLinearCombination<FunctionType> an Eigen scalar type.
 */

namespace Eigen {


template <typename Coefficient, typename FunctionType>
struct NumTraits<GQCP::EvaluableLinearCombination<Coefficient, FunctionType>>:
    public NumTraits<double> {  // permits to get the epsilon, dummy_precision, lowest, highest functions

    using Real = GQCP::EvaluableLinearCombination<Coefficient, FunctionType>;
    using NonInteger = GQCP::EvaluableLinearCombination<Coefficient, FunctionType>;
    using Nested = GQCP::EvaluableLinearCombination<Coefficient, FunctionType>;

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


// Enable the scalar product of a EvaluableLinearCombination<Coefficient, FunctionType> with its own Coefficient (both sides).

template <typename Coefficient, typename FunctionType>
struct ScalarBinaryOpTraits<GQCP::EvaluableLinearCombination<Coefficient, FunctionType>, Coefficient,
                            Eigen::internal::scalar_product_op<GQCP::EvaluableLinearCombination<Coefficient, FunctionType>, Coefficient>> {

    using ReturnType = GQCP::EvaluableLinearCombination<Coefficient, FunctionType>;
};

template <typename Coefficient, typename FunctionType>
struct ScalarBinaryOpTraits<Coefficient, GQCP::EvaluableLinearCombination<Coefficient, FunctionType>,
                            Eigen::internal::scalar_product_op<Coefficient, GQCP::EvaluableLinearCombination<Coefficient, FunctionType>>> {

    using ReturnType = GQCP::EvaluableLinearCombination<Coefficient, FunctionType>;
};


}  // namespace Eigen
