#ifndef GQCG_ATOM_HPP
#define GQCG_ATOM_HPP


#include <iostream>
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
    bool operator==(const GQCG::Atom& other) const;

    /**
     *  @return if this is smaller than @param other, within the @member tolerance_for_comparison for the coordinates
     *
     *  @member atomic_number takes precedence over @member x, over @member y, over @member z
     */
    bool operator<(const GQCG::Atom& other) const;

    /**
     *  Overloading of operator<< for a GQCG::Atom to be used with streams
     */
    friend std::ostream& operator<<(std::ostream& os, const GQCG::Atom& atom);


    // PUBLIC METHODS
    /**
     *  @return if this is equal to @param other, within the given @param tolerance for the coordinates
     */
    bool isEqualTo(const GQCG::Atom& other, double tolerance=Atom::tolerance_for_comparison) const;

    /**
     *  @return if this is smaller than @param other, within the given @param tolerance for the coordinates
     *
     *  @member atomic_number takes precedence over @member x, over @member y, over @member z
     */
    bool isSmallerThan(const GQCG::Atom& other, double tolerance=Atom::tolerance_for_comparison) const;

    /**
     * @return the distance between this and @param other
     */
    double calculateDistance(const GQCG::Atom& other) const;
};


}  // namespace GQCG



#endif  // GQCG_ATOM_HPP
