#ifndef GQCG_ATOM_HPP
#define GQCG_ATOM_HPP


#include <stdlib.h>


namespace GQCG {

/**
 *  A data-holder struct to represent an atom with
 *      - an @member atomic_charge
 *      - coordinates @member x, @member y, @member z
 */
struct Atom {
public:
    size_t atomic_number;
    double x;
    double y;
    double z;

public:
    // CONSTRUCTORS
    /**
     *  Constructor based on a given @param atomic_number and the coordinates @param x, @param y, @param z
     */
    Atom (size_t atomic_number, double x, double y, double z);


    // PUBLIC METHODS
    /**
     *  @return if this is equal to @param other, within the given @param tolerance for the coordinates
     */
    bool isEqualTo(const GQCG::Atom& other, double tolerance=1.0e-08) const;
};


}  // namespace GQCG


#endif  // GQCG_ATOM_HPP
