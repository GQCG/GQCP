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
     *  Constructor based on a given @param atomic_number and the coordinates @param x, @param y, @param z
     */
    Atom (size_t atomic_number, double x, double y, double z);


    // OPERATORS
    /**
     *  @return if this is equal to @param other, within the @member tolerance_for_comparison for the coordinates
     */
    bool operator==(const GQCP::Atom& other) const;

    /**
     *  @return if this is smaller than @param other, within the @member tolerance_for_comparison for the coordinates
     *
     *  @member atomic_number takes precedence over @member x, over @member y, over @member z
     */
    bool operator<(const GQCP::Atom& other) const;

    /**
     *  Overloading of operator<< for a GQCP::Atom to be used with streams
     */
    friend std::ostream& operator<<(std::ostream& os, const GQCP::Atom& atom);


    // PUBLIC METHODS
    /**
     *  @return if this is equal to @param other, within the given @param tolerance for the coordinates
     */
    bool isEqualTo(const GQCP::Atom& other, double tolerance=Atom::tolerance_for_comparison) const;

    /**
     *  @return if this is smaller than @param other, within the given @param tolerance for the coordinates
     *
     *  @member atomic_number takes precedence over @member x, over @member y, over @member z
     */
    bool isSmallerThan(const GQCP::Atom& other, double tolerance=Atom::tolerance_for_comparison) const;

    /**
     * @return the distance between this and @param other
     */
    double calculateDistance(const GQCP::Atom& other) const;
};


}  // namespace GQCP



#endif  // GQCP_ATOM_HPP
