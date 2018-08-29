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
    /**
     *  Constructor based on a given @param atomic_number and the coordinates @param x, @param y, @param z
     */
    Atom (size_t atomic_number, double x, double y, double z);
};


}  // namespace GQCG


#endif  // GQCG_ATOM_HPP
