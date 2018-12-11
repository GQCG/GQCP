#ifndef HoppingMatrix_hpp
#define HoppingMatrix_hpp


#include <Eigen/Dense>


namespace GQCP {


/**
 *  A class that represents a Hubbard hopping matrix
 */
class HoppingMatrix {
private:
    const size_t K;  // the number of lattice sites
    const Eigen::MatrixXd H;  // the Hubbard hopping matrix


public:
    // CONSTRUCTORS
    /**
     *  @param H        the Hubbard hopping matrix
     */
    explicit HoppingMatrix(const Eigen::MatrixXd& H);


    // NAMED CONSTRUCTORS
    /**
     *  @param upper_triangle       the upper triangle (in column-major ordering) that specifies the Hubbard hopping matrix
     *
     *  @return the hopping matrix that corresponds to the given upper triangle
     */
    static HoppingMatrix FromUpperTriangle(const Eigen::VectorXd& upper_triangle);

    /**
     *  @param K        the number of lattice sites
     *
     *  @return a random hopping matrix with elements distributed uniformly in [-1.0, 1.0]
     */
    static HoppingMatrix Random(size_t K);


    // PUBLIC METHODS
    /**
     *  @return the Hubbard hopping matrix
     */
    const Eigen::MatrixXd& asMatrix() const { return this->H; }

    /**
     *  @return the number of lattice sites corresponding to the Hubbard hopping matrix
     */
    size_t numberOfLatticeSites() const { return this->K; }
};


}  // namespace GQCP


#endif /* HoppingMatrix_hpp */
