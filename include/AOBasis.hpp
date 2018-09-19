#ifndef GQCG_AOBASIS_HPP
#define GQCG_AOBASIS_HPP


#include "Atom.hpp"
#include "Molecule.hpp"

#include <Eigen/Dense>
#include <libint2.hpp>




namespace GQCG {


/**
 *  A class that represents an atomic orbital basis
 */
class AOBasis {
private:
    const std::vector<GQCG::Atom> atoms;
    const libint2::BasisSet basis_functions;
    const size_t number_of_basis_functions;

public:
    // CONSTRUCTOR
    AOBasis(const GQCG::Molecule& molecule, std::string basis_set);


    // GETTERS
    size_t get_number_of_basis_functions() const { return this->number_of_basis_functions; }


    // FRIEND CLASSES
    friend class LibintCommunicator;
};


}  // namespace GQCG


#endif  // GQCG_AOBASIS_HPP
