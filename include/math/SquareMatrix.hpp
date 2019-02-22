#ifndef SquareMatrix_hpp
#define SquareMatrix_hpp


#include "typedefs.hpp"


namespace GQCP {


/**
 *  A square extension of the Matrix class
 *
 *  @tparam Scalar      the scalar type
 */
template<typename Scalar>
class SquareMatrix : public Matrix<Scalar> {
public:

    /*
     *  CONSTRUCTORS
     */

    /**
     *  Default constructor
     */
    SquareMatrix() : Matrix<Scalar>() {}


    /**
     *  A basic constructor that checks if the given matrix is square
     *
     *  @param matrix       the matrix that should be square
     */
    SquareMatrix(const Matrix<Scalar>& matrix) :
        Matrix<Scalar>(matrix)  // the compiler should call the move constructor here
    {
        // Check if the given matrix is square
        if (this->cols() != this->rows()) {
            throw std::invalid_argument("The given matrix is not square.");
        }
    }


    /*
     *  GETTERS
     */

    size_t get_dim() const {
        return this->cols();  // equals this->rows()
    }
};


}  // namespace GQCP


#endif /* SquareMatrix_hpp */
