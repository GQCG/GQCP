// This file is part of GQCG-gqcp.
// 
// Copyright (C) 2017-2019  the GQCG developers
// 
// GQCG-gqcp is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// GQCG-gqcp is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-gqcp.  If not, see <http://www.gnu.org/licenses/>.
// 
#ifndef GQCP_ATOM_HPP
#define GQCP_ATOM_HPP


#include "Mathematical/Matrix.hpp"

#include <iostream>
#include <stdlib.h>


namespace GQCP {


/**
 *  A class that represents a nucleus: it has a charge and a position in space
 */
class Nucleus {
public:
    size_t atomic_number;
    Vector<double, 3> position;  // in bohr

    static constexpr double tolerance_for_comparison = 1.0e-08;  // in bohr

public:
    // CONSTRUCTORS
    /**
     *  @param atomic_number        the atomic number (Z) of the nucleus
     *  @param x                    the x-position of the nucleus in bohr
     *  @param y                    the y-position of the nucleus in bohr
     *  @param z                    the z-position of the nucleus in bohr
     */
    Nucleus(size_t atomic_number, double x, double y, double z);

    /**
     *  Default constructor, creating a 'ghost' nucleus (i.e. Bq) in the origin
     */
    Nucleus();


    // OPERATORS
    /**
     *  @param other        the other nucleus
     *
     *  @return if this nucleus is equal to the other, within a default tolerance for the coordinates
     */
    bool operator==(const Nucleus& other) const;

    /**
     *  @param other        the other nucleus
     *
     *  @return if this nucleus is not equal to the other, within a default tolerance for the coordinates
     */
    bool operator!=(const Nucleus& other) const;

    /**
     *  A custom implementation for the comparison (and thus ordening) of nuclei. The atomic_number takes precedence over the x-coordinate, which takes precedence over the y-coordinate, which in turn takes precedence over the z-coordinate
     *
     *  @param other        the other nucleus
     *
     *  @return if this nucleus is 'smaller' than the other, within a default tolerance for the coordinates
     */
    bool operator<(const Nucleus& other) const;

    /**
     *  Overloading of operator<< for a Nucleus to be used with ostreams
     *
     *  @param os       the output stream to which the nucleus should be concatenated
     *  @param nucleus     the nucleus which should be concatenated to the output stream
     *
     *  @return the updated output stream
     */
    friend std::ostream& operator<<(std::ostream& os, const Nucleus& nucleus);


    // PUBLIC METHODS
    /**
     *  @param other        the other nucleus
     *  @param tolerance    the tolerance for equality of positions
     *
     *  @return if this nucleus is equal to the other
     */
    bool isEqualTo(const Nucleus& other, double tolerance=Nucleus::tolerance_for_comparison) const;

    /**
     *  A custom implementation for the comparison (and thus ordening) of nuclei. The atomic_number takes precedence over the x-coordinate, which takes precedence over the y-coordinate, which in turn takes precedence over the z-coordinate
     *
     *  @param other        the other nucleus
     *  @param tolerance    the tolerance for equality of positions
     *
     *  @return if this nucleus is 'smaller' than the other, within a default tolerance for the coordinates
     */
    bool isSmallerThan(const Nucleus& other, double tolerance=Nucleus::tolerance_for_comparison) const;

    /**
     *  @param other        the other nucleus
     *
     *  @return the Euclidian distance between this nucleus and the other
     */
    double calculateDistance(const Nucleus& other) const;
};


}  // namespace GQCP



#endif  // GQCP_ATOM_HPP
