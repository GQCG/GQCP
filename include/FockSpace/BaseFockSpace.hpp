#ifndef GQCG_BASEFOCKSPACE_HPP
#define GQCG_BASEFOCKSPACE_HPP


#include "ONV.hpp"
#include "FockSpace/FockSpaceType.hpp"
#include "common.hpp"



namespace GQCG {


/**
 *  A base class for the Fock space
 *  Interfacing requires the Fock space to generate an ONV from a given address
 *  transform a given ONV into the next ONV (in the full or selected space)
 *  and retrieve the address of a given ONV in the space
 */
class BaseFockSpace {
protected:
    const size_t K;  // number of spatial orbitals
    const size_t dim;  // dimension of the Fock space


    // PROTECTED CONSTRUCTORS
    /**
     *  Protected constructor given a @param K and @param dim
     */
    explicit BaseFockSpace(size_t K, size_t dim);


public:
    // DESTRUCTOR
    /**
     *  Provide a pure virtual destructor to make the class abstract
     */
    virtual ~BaseFockSpace() = 0;


    // GETTERS
    size_t get_dimension() const { return dim; }
    size_t get_K() const { return K; }
    virtual FockSpaceType get_type() const = 0;

    // PUBLIC METHODS
    /**
     *  Creates a Hartree-Fock coefficient expansion (single Slater expansion of the first configuration in the Fock space)
     */
    Eigen::VectorXd HartreeFockExpansion();
};


}  // namespace GQCG


#endif  // GQCG_BASEFOCKSPACE_HPP
