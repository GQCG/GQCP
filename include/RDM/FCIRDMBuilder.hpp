#ifndef GQCG_FCIRDMBUILDER_HPP
#define GQCG_FCIRDMBUILDER_HPP


#include "FockSpace/FockSpaceProduct.hpp"
#include "RDM/RDMBuilder.hpp"
#include "RDM/RDMs.hpp"


namespace GQCG {


/**
 *  FCIRDMBuilder is a class for the calculation of a density matrix from a given wave function
 *  or coefficient expansion in a the Full ci Fock space
 */
class FCIRDMBuilder : public RDMBuilder {
    FockSpaceProduct fock_space;  // fock space containing the alpha and beta Fock space


public:
    // CONSTRUCTOR
    explicit FCIRDMBuilder(FockSpaceProduct fock_space);


    // DESTRUCTOR
    ~FCIRDMBuilder() = default;


    // OVERRIDDEN PUBLIC METHODS
    /**
     *  @return all 1-RDMs from a coefficient vector @param x
     */
    OneRDMs calculate1RDMs(const Eigen::VectorXd &x) override;

    /**
     *  @return all 2-RDMs from a coefficient vector @param x
     */
    TwoRDMs calculate2RDMs(const Eigen::VectorXd &x) override;

    /**
     *  @return the Fock space of the RDMBuilder
     */
    BaseFockSpace* get_fock_space() override { return &fock_space; }
};


}  // namespace GQCG


#endif  // GQCG_FCIRDMBUILDER_HPP
