#ifndef BasisFunction_hpp
#define BasisFunction_hpp


#include "Basis/CartesianGTO.hpp"
#include "math/LinearCombination.hpp"



namespace GQCP {


/**
 *  A class that represents a basis function and can be evaluated at a point in Euclidean space
 */
class BasisFunction : public LinearCombination<double, CartesianGTO> {
private:
    double N_total;  // the total normalization factor


public:
    using Base = LinearCombination<double, CartesianGTO>;


public:
    // CONSTRUCTORS
    /**
     *  @param lc       a linear combination of CartesianGTOs
     */
    BasisFunction(const Base& lc);


    // PUBLIC METHODS
    /**
     *  @return the total normalization factor of this basis function
     */
    double calculateNormalizationFactor() const;  // TODO: implement
};


}  // namespace GQCP


#endif  /* BasisFunction_hpp */
