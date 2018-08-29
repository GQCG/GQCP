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
     *  A standard constructor to be able to use emplace_back
     * @param atomic_number
     * @param x
     * @param y
     * @param z
     */
    Atom (size_t atomic_number, double x, double y, double z) :
        atomic_number (atomic_number),
        x (x),
        y (y),
        z (z)
    {}
};


}  // namespace GQCG


#endif  // GQCG_ATOM_HPP
