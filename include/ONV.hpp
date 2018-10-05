#ifndef GQCG_ONV_HPP
#define GQCG_ONV_HPP


#include "common.hpp"



namespace GQCG {


class ONV {
private:
    const size_t K;  // number of spatial orbitals
    const size_t N;  // number of electrons
    size_t unsigned_representation;  // unsigned representation
    size_t_sptr occupation_indexes;  // the occupied orbital electron indexes


    // PRIVATE METHODS
    /**
     *  Extracts the positions of the set bits from the representation and places them in an array
     */
    void update();

    /**
     *  Tests whether the assigned amount of electrons N and the amount of set bits in the representation match
     */
    void is_compatible(size_t l);



public:
    // CONSTRUCTORS
    /**
     *  Constructor from a @param K orbitals, N electrons and a representation for the ONV
     */
    ONV(size_t K, size_t N, size_t unsigned_representation);


    // OPERATORS
    /**
     *  Overloading of operator<< for a GQCG::ONV to be used with streams
     */
    friend std::ostream& operator<<(std::ostream& os, const GQCG::ONV& onv);

    /**
     *  @return if this equals @param other
     */
    bool operator==(ONV& other) const {
        return this->unsigned_representation == other.unsigned_representation && this->K == other.K;  // this ensures that N, K and representation are equal
    }

    /**
     *  @return if this is not equal to @param other
     */
    bool operator!=(ONV& other) const {
        return !(this->operator==(other));
    }

    // GETTERS & SETTERS
    void set_representation(size_t unsigned_representation);
    size_t get_urepresentation(){ return unsigned_representation;}

    /**
     *  @return occupied orbital based on the electron index
     */
    size_t get_occupied_orbital(size_t electron_index);


    // PUBLIC METHODS (PREVIOUS SPIN STRING FUNCTIONALITY)
    /**
     *  @return if the index is occupied (i.e. 1) for the @param p-th spatial orbital, starting from 0
     *  @param p is the lexical index (i.e. read from right to left)
     */
    bool isOccupied(size_t p) const;

    /**
     *  @return if we can apply the annihilation operator (i.e. 1->0) for the @param p-th spatial orbital, starting from 0
     *  performs an annihilation @param p
     *
     *  !!! IMPORTANT: does not update the occupation array if required call "update()" !!!
     *  !!! IMPORTANT: performs the annihilation in place !!!
     */
    bool annihilate(size_t p);

    /**
     *  @return if we can apply the annihilation operator (i.e. 1->0) for the @param p-th spatial orbital, starting from 0
     *  @param p is the lexical index (i.e. read from right to left)
     *
     *  Furthermore, the @param sign is changed according to the sign change (+1 or -1) of the spin string after annihilation.
     *
     *  !!! IMPORTANT: does not update the occupation array if required call "update()" !!!
     *  !!! IMPORTANT: performs the annihilation in place !!!
     */
    bool annihilate(size_t p, int& sign);

    /**
     * @return if we can apply the creation operator (i.e. 0->1) for the @param p-th spatial orbital, starting from 0
     *  performs a creation @param p
     *
     *  !!! IMPORTANT: does not update the occupation array if required call "update()" !!!
     *  !!! IMPORTANT: performs the creation in place !!!
     */
    bool create(size_t p);

    /**
     *  @return if we can apply the creation operator (i.e. 0->1) for the @param p-th spatial orbital, starting from 0
     *  @param p is the lexical index (i.e. read from right to left)
     *
     *  Furthermore, the @param sign is changed according to the sign change (+1 or -1) of the spin string after annihilation.
     *
     *  !!! IMPORTANT: does not update the occupation array if required call "update()" !!!
     *  !!! IMPORTANT: performs the creation in place !!!
     */
    bool create(size_t p, int& sign);

    /**
     *  @return the phase factor (+1 or -1) that arises by applying an annihilation or creation operator on orbital @param p, starting from 0 in reverse lexical ordering.
     *
     *  Let's say that there are m electrons in the orbitals up to p (not included). If m is even, the phase factor is (+1) and if m is odd, the phase factor is (-1), since electrons are fermions.
     */
    int operatorPhaseFactor(size_t p) const;


    /**
     *  @return the representation of a slice (i.e. a subset) of the spin string between @param index_start (included)
     *  and @param index_end (not included).
     *
     *  Both @param index_start and @param index_end are 'lexical' (i.e. from right to left), which means that the slice
     *  occurs 'lexically' as well (i.e. from right to left).
     *
     *      Example:
     *          "010011".slice(1, 4) => "01[001]1" -> "001"
     *
     */
    size_t slice(size_t index_start, size_t index_end) const;


    // FRIEND CLASSES
    friend class FockSpace;
};


}  // namespace GQCG

#endif //GQCG_ONV_HPP
