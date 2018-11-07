// This file is part of GQCG-gqcp.
// 
// Copyright (C) 2017-2018  the GQCG developers
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


#include <iostream>
#include <stdlib.h>

#include <Eigen/Dense>


namespace GQCP {

/**
 *  A data-holder struct to represent an atom with
 *      - an @member atomic_charge
 *      - coordinates @member x, @member y, @member z
 */
struct Atom {
public:
    size_t atomic_number;
    Eigen::Vector3d position;

    static constexpr double tolerance_for_comparison = 1.0e-08;

public:
    // CONSTRUCTORS
    /**
     *  @param atomic_number        the atomic number (Z) of the atom
     *  @param x                    the x-position of the atom
     *  @param y                    the y-position of the atom
     *  @param z                    the z-position of the atom
     */
    Atom(size_t atomic_number, double x, double y, double z);


    // OPERATORS
    /**
     *  @param other        the other atom
     *
     *  @return if this atom is equal to the other, within a default tolerance for the coordinates
     */
    bool operator==(const GQCP::Atom& other) const;

    /**
     *  A custom implementation for the comparison (and thus ordening) of atoms. The atomic_number takes precedence over the x-coordinate, which takes precedence over the y-coordinate, which in turn takes precedence over the z-coordinate
     *
     *  @param other        the other atom
     *
     *  @return if this atom is 'smaller' than the other, within a default tolerance for the coordinates
     */
    bool operator<(const GQCP::Atom& other) const;

    /**
     *  Overloading of operator<< for a GQCP::Atom to be used with ostreams
     *
     *  @param os       the output stream to which the atom should be concatenated
     *  @param atom     the atom which should be concatenated to the output stream
     *
     *  @return the updated output stream
     */
    friend std::ostream& operator<<(std::ostream& os, const GQCP::Atom& atom);


    // PUBLIC METHODS
    /**
     *  @param other        the other atom
     *  @param tolerance    the tolerance for equality of positions
     *
     *  @return if this atom is equal to the other
     */
    bool isEqualTo(const GQCP::Atom& other, double tolerance=Atom::tolerance_for_comparison) const;

    /**
     *  A custom implementation for the comparison (and thus ordening) of atoms. The atomic_number takes precedence over the x-coordinate, which takes precedence over the y-coordinate, which in turn takes precedence over the z-coordinate
     *
     *  @param other        the other atom
     *  @param tolerance    the tolerance for equality of positions
     *
     *  @return if this atom is 'smaller' than the other, within a default tolerance for the coordinates
     */
    bool isSmallerThan(const GQCP::Atom& other, double tolerance=Atom::tolerance_for_comparison) const;

    /**
     *  @param other        the other atom
     *
     *  @return the Euclidian distance between this atom and the other
     */
    double calculateDistance(const GQCP::Atom& other) const;
};


}  // namespace GQCP



#endif  // GQCP_ATOM_HPP
