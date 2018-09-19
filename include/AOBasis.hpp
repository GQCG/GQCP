#ifndef GQCG_AOBASIS_HPP
#define GQCG_AOBASIS_HPP


#include "Atom.hpp"
#include "Molecule.hpp"

#include <Eigen/Dense>
#include <libint2.hpp>




namespace GQCG {


/**
 *  A class that represents an atomic orbital basis, with the members
 *      - @member basis_functions
 *      - @member S, which is the overlap matrix of the @member basis_functions
 */
class AOBasis {
private:
    const std::vector<GQCG::Atom> atoms;
    const libint2::BasisSet basis_functions;

public:
    // CONSTRUCTOR
    AOBasis(const GQCG::Molecule& molecule, std::string basis_set);


    // FRIEND CLASSES
    friend class LibintCommunicator;
};


}  // namespace GQCG


#endif  // GQCG_AOBASIS_HPP
