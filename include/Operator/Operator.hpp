#ifndef IntegralParameters_hpp
#define IntegralParameters_hpp


#include "typedefs.hpp"
#include "math/SquareMatrix.hpp"


namespace GQCP {

/**
 *  An interface for second-quantized operators: they should implement the transformation formulas for their matrix representations in an orbital basis
 */
template<typename Scalar>
class Operator {
public:

    /*
     *  PUBLIC PURE VIRTUAL METHODS
     */

    /**
     *  In-place transform the matrix representation of the operator
     *
     *  @param T    the transformation matrix between the old and the new orbital basis, it is used as
     *      b' = b T ,
     *   in which the basis functions are collected as elements of a row vector b
     */
    virtual void transform(const SquareMatrix<Scalar>& T) = 0;


    /*
     *  PUBLIC METHODS
     */

    /**
     *  In-place rotate the matrix representation of the operator
     *
     *  @param U     the unitary transformation (i.e. rotation) matrix, see transform() for how the transformation matrix between the two bases should be represented
     */
    void rotate(const SquareMatrix<Scalar>& U) {

        // Check if the given matrix is actually unitary
        if (!U.isUnitary(1.0e-12)) {
            throw std::invalid_argument("The given transformation matrix is not unitary.");
        }

        this->transform(U);
    }


    // Normally, I would have liked to put in the 'only-if-double' Jacobi rotation formula, but a function can't be and templated and virtual
};


}


#endif /* IntegralParameters_hpp */
