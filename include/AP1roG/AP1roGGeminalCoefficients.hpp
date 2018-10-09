#ifndef AP1roGGeminalCoefficients_hpp
#define AP1roGGeminalCoefficients_hpp


#include <Eigen/Dense>


namespace GQCG {


class AP1roGGeminalCoefficients {
private:
    size_t N_P;  // the number of electron pairs (= the number of geminals)
    size_t K;  // the number of orbitals
    Eigen::VectorXd g;


public:
    // CONSTRUCTORS
    /**
     *  Default constructor setting everything to zero
     */
    AP1roGGeminalCoefficients();
    
    /**
     *  Constructor setting the geminal coefficients to zero, based on the number of orbitals @param K and number of electron pairs @param N_P
     */
    AP1roGGeminalCoefficients(size_t N_P, size_t K);

    /**
     *  Constructor based on given geminal coefficients @param g, which are in vector (row-major) form
     */
    AP1roGGeminalCoefficients(const Eigen::VectorXd& g, size_t N_P, size_t K);


    // METHODS
    /**
     *  @return the geminal coefficients in vector form
     */
    Eigen::VectorXd asVector() const { return this->g; }

    /**
     *  Construct and @return the geminal coefficients in matrix form
     */
    Eigen::MatrixXd asMatrix() const;
};


}  // namespace GQCG



#endif /* AP1roGGeminalCoefficients_hpp */
