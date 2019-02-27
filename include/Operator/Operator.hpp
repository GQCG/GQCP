#ifndef Operator_hpp
#define Operator_hpp


#include "typedefs.hpp"
#include "math/SquareMatrix.hpp"


namespace GQCP {

/**
 *  An interface for second-quantized operators: they should implement the transformation formulas for their matrix representations in an orbital basis
 *
 *  Since Operator::rotate() is implemented in the base class using a derived-class transform(), we use CRTP as static polymorphism
 */
template <typename DerivedOperator>
class Operator {
public:

    /**
     *  @return this as a DerivedOperator (done at compile time)
     */
    DerivedOperator& derived() { return static_cast<DerivedOperator&>(*this); }

    /**
     *  @return this as a const DerivedOperator (done at compile time)
     */
    const DerivedOperator& derived() const { return static_cast<DerivedOperator&>(*this); }


    /*
     *  PUBLIC METHODS
     */

    /**
     *  In-place rotate the matrix representation of the operator
     *
     *  @param U     the unitary transformation (i.e. rotation) matrix, see transform() for how the transformation matrix between the two bases should be represented
     */
    template<typename Scalar>
    void rotate(const SquareMatrix<Scalar>& U) {

        // Check if the given matrix is actually unitary
        if (!U.isUnitary(1.0e-12)) {
            throw std::invalid_argument("The given transformation matrix is not unitary.");
        }

        this->derived().transform(U);
    }


    // Normally, I would have liked to put in the 'only-if-double' Jacobi rotation formula, but a function can't be and templated and virtual
};


}  // namespace GQCP


#endif  /* Operator_hpp */
