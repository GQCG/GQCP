#ifndef AP1roGVariables_hpp
#define AP1roGVariables_hpp


#include "geminals/BaseAPIGVariables.hpp"


namespace GQCP {


class AP1roGVariables : public BaseAPIGVariables {
public:
    // CONSTRUCTORS
    /**
     *  Default constructor setting everything to zero
     */
    AP1roGVariables();  // default constructor needed

    /**
     *  @param x        the variables in a vector representation that is in row-major storage
     *
     *  @param N_P      the number of electron pairs (= the number of geminals)
     *  @param K        the number of spatial orbitals
     */
    AP1roGVariables(const Eigen::VectorXd& x, size_t N_P, size_t K);

    /**
     *  Constructor that sets the variables to zero
     *
     *  @param N_P      the number of electron pairs (= the number of geminals)
     *  @param K        the number of spatial orbitals
     */
    AP1roGVariables(size_t N_P, size_t K);


    // STATIC PUBLIC METHODS
    /**
     *  @param N_P      the number of electron pairs (= the number of geminals)
     *  @param K        the number of spatial orbitals
     *
     *  @return the total number of variables
     */
    static size_t numberOfVariables(size_t N_P, size_t K);


    // PUBLIC METHODS
    /**
     *  @return the variables in matrix form
     */
    virtual Eigen::MatrixXd asMatrix() const override;

    /**
     *  @param vector_index     the vector index of the variable
     *
     *  @return the major (geminal, non-contiguous) index i (i.e. the subscript) in the matrix of the variables. Note that i is in [0 ... N_P[
     */
    size_t matrixIndexMajor(size_t vector_index) const override;

    /**
     *  @param vector_index     the vector index of the variable
     *
     *  @return the minor (virtual orbital, contiguous) index a (i.e. the subscript) in the matrix of the variables. Note that a is in [N_P ... K[
     */
    size_t matrixIndexMinor(size_t vector_index) const override;

    /**
     *  @param i        the major (geminal) index (changes in i are not contiguous)
     *  @param a        the minor (virtual orbital) index (changes in a are contiguous)
     *
     *  @return the vector index of the variable X_i^a
     */
    size_t vectorIndex(size_t i, size_t a) const override;
};



}  // namespace GQCP




#endif /* AP1roGVariables_hpp */
