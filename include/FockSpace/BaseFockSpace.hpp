#ifndef GQCG_BASEFOCKSPACE_HPP
#define GQCG_BASEFOCKSPACE_HPP


#include "ONV.hpp"
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
    size_t dim;  // dimension of the Fock space


    // PROTECTED CONSTRUCTORS
    /**
     *  Protected constructor given a @param K
     */
    explicit BaseFockSpace(size_t K);


public:
    // DESTRUCTOR
    /**
     *  Provide a pure virtual destructor to make the class abstract
     */
    virtual ~BaseFockSpace() = 0;


    // GETTERS
    size_t get_dimension() const { return dim; }
    size_t get_K() const { return K; }


    // PUBLIC METHODS
    /**
     *  Creates a Hartree-Fock guess
     */
    Eigen::VectorXd hartreeFockGuess();

    // PURE VIRTUAL PUBLIC METHODS
    /**
     *  @return the ONV with the corresponding @param address in the considered space
     */
    virtual ONV get_ONV(size_t address) = 0;

    /**
     *  sets @param ONV to the next ONV in the space
     */
    virtual void setNext(ONV& onv) = 0;

    /**
     *  @return the Fock space address (i.e. the ordering number) of the @param onv in reverse lexical ordering, in the fock space.
     */
    virtual size_t getAddress(ONV& onv) = 0;
};


}  // namespace GQCG


#endif  // GQCG_BASEFOCKSPACE_HPP
