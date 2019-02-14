#ifndef SquareFourIndexTensor_hpp
#define SquareFourIndexTensor_hpp


#include "typedefs.hpp"



namespace GQCP {


/**
 *  A square extension of a 4-index tensor class
 *
 *  @tparam Scalar      the scalar type
 */
template<typename Scalar>
class SquareFourIndexTensor: public FourIndexTensor<Scalar> {
public:

    // CONSTRUCTORS
    /**
     *  Default constructor
     */
    SquareFourIndexTensor() : FourIndexTensor<Scalar>() {}


    /**
     *  A basic constructor that checks if the given tensor is 'square'
     *
     *  @param tensor       the tensor that should be 'square'
     */
    SquareFourIndexTensor(const FourIndexTensor<Scalar>& tensor) :
        FourIndexTensor<Scalar>(tensor)
    {
        // Check if the given tensor is 'square'
        auto dims = tensor.dimensions();
        if ((dims[0] != dims[1]) || (dims[1] != dims[2]) || (dims[2] != dims[3]) ) {
            throw std::invalid_argument("The given tensor should have equal dimensions in every rank.");
        }
    }


    // GETTERS
    size_t get_dim() const {
        return this->dimension(0);  // all tensor dimensions are equal because of the constructor
    }
};


}  // namespace GQCP


#endif /* FourIndexTensor_hpp */
