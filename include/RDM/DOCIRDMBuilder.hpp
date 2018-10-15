#ifndef GQCG_DOCIRDMBUILDER_HPP
#define GQCG_DOCIRDMBUILDER_HPP


namespace GQCG {


/**
 *  DOCIRDMBuilder is a class for the calculation of a density matrix from a given wave function
 *  or coefficient expansion in a doubly occupied or single Fock space
 */
class DOCIRDMBuilder {
    FockSpace fock_space;  // both the alpha and beta Fock space


public:
    // CONSTRUCTOR
    explicit DOCIRDMBuilder(FockSpace fock_space);


    // DESTRUCTOR
    ~DOCIRDMBuilder() = default;


    // PURE VIRTUAL PUBLIC METHODS
    /**
     *  @return 1RDM from a coefficient vector @param x
     */
    Eigen::MatrixXd construct1RDM(const Eigen::VectorXd& x) override;

    /**
     *  @return 2RDM from a coefficient vector @param x
     */
    Eigen::Tensor<double, 4> construct2RDM(const Eigen::VectorXd& x) override;

    /**
     *  @return the Fock space of the RDMBuilder
     */
    virtual BaseFockSpace* get_fock_space() override { return &fock_space; }
};


}  // namespace GQCG


#endif //GQCG_DOCIRDMBUILDER_HPP
