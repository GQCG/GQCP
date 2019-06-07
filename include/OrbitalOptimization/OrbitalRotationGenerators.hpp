#ifndef GQCP_ORBITALROTATIONGENERATORS_HPP
#define GQCP_ORBITALROTATIONGENERATORS_HPP


#include "math/SquareMatrix.hpp"


namespace GQCP {


/**
 *  A class that represents the orbital rotation generators kappa
 */
class OrbitalRotationGenerators {
private:
    size_t number_of_spatial_orbitals;  // the number of spatial orbitals that can be rotated using these orbital rotation generators
    VectorX<double> kappa_vector;  // the strict upper/lower triangle of the kappa matrix


public:
    // CONSTRUCTORS

    /**
     *  @param  kappa_vector        the orbital rotation generators represented as a vector that corresponds to the strict upper/lower triangle of the kappa matrix
     */
    OrbitalRotationGenerators(const VectorX<double>& kappa_vector);

    /**
     *  @param  kappa_vector        the orbital rotation generators represented as the full antisymmetric matrix kappa
     */
    OrbitalRotationGenerators(const SquareMatrix<double>& kappa_matrix);


    // NAMED CONSTRUCTORS

    /**
     *  Construct orbital rotation generators by adding redundant (i.e. 0) occupied-virtual and virtual-virtual generators to the given occupied-occupied generators
     * 
     *  @param o_o_generators       the occupied-occupied orbital rotation generators
     *  @param K                    the total number of spatial orbitals
     * 
     *  @return 'full' orbital rotation generators from the given occupied-occupied generators
     */
    static OrbitalRotationGenerators FromOccOcc(const OrbitalRotationGenerators& o_o_generators, const size_t K);


    // PUBLIC METHODS

    /**
     *  @return the orbital rotation generators as the strict upper/lower triangle of the kappa matrix
     */
    const VectorX<double>& asVector() const { return this->kappa_vector; }

    /**
     *  @return the antisymmetric orbital rotation generator matrix kappa
     */
    SquareMatrix<double> asMatrix() const;

    /**
     *  @return the unitary matrix that corresponds to these orbital rotation generators, i.e. exp(-kappa)
     */
    SquareMatrix<double> calculateRotationMatrix() const;

    /*
     *  @return the number of spatial orbitals that can be rotated using these orbital rotation generators
     */
    size_t numberOfSpatialOrbitals() const { return this->number_of_spatial_orbitals; }
};


}  // namespace GQCP



#endif  /* GQCP_ORBITALROTATIONGENERATORS_HPP */
