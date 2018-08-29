#include "Atom.hpp"


namespace GQCG {


/*
 *  CONSTRUCTORS
 */

/**
 *  Constructor based on a given @param atomic_number and the coordinates @param x, @param y, @param z
 */
Atom::Atom(size_t atomic_number, double x, double y, double z) :
    atomic_number (atomic_number),
    x (x),
    y (y),
    z (z)
{}



}  // namespace GQCG
