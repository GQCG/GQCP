#ifndef SquareMatrix_hpp
#define SquareMatrix_hpp


#include "typedefs.hpp"


namespace GQCP {


/**
 *  A class that is a square extension of the MatrixX class
 *
 *  @tparam Scalar      the scalar type
 */
template<typename Scalar>
class SquareMatrixX : public MatrixX<Scalar> {
public:

    // CONSTRUCTORS
    /**
     *  (Move) constructor from a matrix representation
     *
     *  @param matrix       the matrix representation
     */
    SquareMatrixX(const MatrixX<Scalar>& matrix) :
        MatrixX<Scalar>(matrix)
    {
        if (this->cols() != this->rows()) {
            throw std::invalid_argument("The given matrix is not square.");
        }
    }
};


}  // namespace GQCP


#endif /* SquareMatrix_hpp */
