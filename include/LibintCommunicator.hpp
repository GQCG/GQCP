#ifndef GQCG_LIBINTCOMMUNICATOR_HPP
#define GQCG_LIBINTCOMMUNICATOR_HPP


#include "AOBasis.hpp"
#include "Molecule.hpp"
#include "Operator/OneElectronOperator.hpp"
#include "Operator/TwoElectronOperator.hpp"

#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>
#include <libint2.hpp>


namespace GQCG {


/**
 *  A singleton class that takes care of interfacing with the Libint2 (version >2.2.0) C++ API
 *
 *  Singleton class template from (https://stackoverflow.com/a/1008289)
 */
class LibintCommunicator {
private:
    /**
     *  Private constructor as required by the singleton class design
     */
    LibintCommunicator();

    /**
     *  Private destructor as required by the singleton class design
     */
    ~LibintCommunicator();

public:
    /**
     *  @return the static singleton instance
     */
    static LibintCommunicator& get();


    /**
     *  Remove the public copy constructor and a public assignment operator
     */
    LibintCommunicator(LibintCommunicator const& libint_communicator) = delete;
    void operator=(LibintCommunicator const& libint_communicator) = delete;


    // PUBLIC METHODS
    /**
     *  @return a std::vector<libint2::Atom> based on a given std::vector<GQCG::Atom> @param atoms
     */
    std::vector<libint2::Atom> interface(const std::vector<GQCG::Atom>& atoms) const;


    /**
     *  @return the OneElectronOperator corresponding to the matrix representation of @param operator_type in the given
     *  @param ao_basis. The corresponding @param molecule is also given as an argument, to be able to access
     */
    GQCG::OneElectronOperator calculateOneElectronIntegrals(libint2::Operator operator_type, const GQCG::AOBasis& ao_basis) const;

    /**
     *  @return the TwoElectronOperator corresponding to the matrix representation of @param operator_type in the given
     *  @param ao_basis
     */
    GQCG::TwoElectronOperator calculateTwoElectronIntegrals(libint2::Operator operator_type, const GQCG::AOBasis& ao_basis) const;
};


}  // namespace GQCG


#endif  // GQCG_LIBINTCOMMUNICATOR_HPP
