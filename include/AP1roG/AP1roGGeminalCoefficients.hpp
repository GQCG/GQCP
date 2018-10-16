#ifndef AP1roGGeminalCoefficients_hpp
#define AP1roGGeminalCoefficients_hpp


#include <Eigen/Dense>


namespace GQCG {


class AP1roGGeminalCoefficients {
private:
    size_t N_P;  // the number of electron pairs (= the number of geminals)
    size_t K;  // the number of orbitals
    Eigen::VectorXd g;  // the geminal coefficients stored in a row-major form


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


    // OPERATORS
    /**
     *  @return the geminal coefficient g_mu, in which @param mu is a 'vector index'
     */
    double operator()(size_t mu) const;

    /**
     *  @ return the geminal coefficient G_i^a, in which @param i is the major index (changes in i are not contiguous) and @param a is the minor index (changes in a are contiguous)
     */
    double operator()(size_t i, size_t a) const;


    // GETTERS
    size_t get_N_P() const { return this->N_P; }
    size_t get_K() const { return this->K; }


    // METHODS
    /**
     *  @return the geminal coefficients in vector form
     */
    Eigen::VectorXd asVector() const { return this->g; }

    /**
     *  Construct and @return the geminal coefficients in matrix form
     */
    Eigen::MatrixXd asMatrix() const;

    /**
     *  @return the number of free geminal coefficients given a number of geminals @param N_P and a number of spatial orbitals @param K
     */
    static size_t numberOfGeminalCoefficients(size_t N_P, size_t K);

    /**
     *  For a geminal coefficient g_mu, return its major index in the matrix of geminal coefficients.
     *
     *      Note that:
     *          - the major index is i (i.e. the subscript), since changes in i are not contiguous
     *          - i is in [0 ... N_P[
     */
    size_t matrixIndexMajor(size_t vector_index) const;

    /**
     *  For a geminal coefficient g_mu, return its minor index in the matrix of geminal coefficients.
     *
     *      Note that:
     *          - the minor index is a (i.e. the superscript), since changes in a are contiguous
     *          - a is in [N_P ... K[
     */
    size_t matrixIndexMinor(size_t vector_index) const;

    /**
     *  For a geminal coefficient G_i^a, return its index in the vector of geminal coefficients.
     *
     *      Note that
     *          - i is in [0 ... N_P[       is the 'major' index (i.e. changes in i are not contiguous)
     *          - a is in [N_P ... K[       is the 'minor' index (i.e. changes in a are contiguous)
     */
    size_t vectorIndex(size_t i, size_t a) const;
};


}  // namespace GQCG



#endif /* AP1roGGeminalCoefficients_hpp */
