#include <bitset>
#include "ONV.hpp"


namespace GQCG {


/*
 *  PRIVATE METHODS
 */

/**
 *  Extracts the positions of the set bits from the representation and places them in an array
 */
void ONV::update(){
    size_t l = this->unsigned_representation;
    int i = 0;
    while(l != 0){
        this->occupation_indexes(i) = __builtin_ctzl(l);
        i++;
        l ^= (l & -l);
    }
}

/**
 *  Tests whether the assigned amount of electrons N and the amount of set bits in the representation match
 */
void ONV::is_compatible(size_t l){
    int i = 0;
    while(l != 0){
        if(i == this->N){
            throw std::invalid_argument("representation and electron count are not compatible");
        }
        i++;
        l ^= (l & -l);
    }
}



/*
 *  CONSTRUCTORS
 */

/**
 *  Constructor from a @param K orbitals, N electrons and a representation for the ONV
 */
ONV::ONV(size_t K, size_t N, size_t representation): K(K), N(N), unsigned_representation(representation){
    occupation_indexes = Eigen::VectorXd::Zero(N);
    is_compatible(this->unsigned_representation);  // throws error if the constructor parameters are not compatible;
    update();
}



/*
 *  OPERATORS
 */

/**
 *  Overloading of operator<< for a GQCG::ONV to be used with streams
 */
std::ostream& operator<<(std::ostream& os, const GQCG::ONV& onv) {
    const char* begin = reinterpret_cast<const char*>(&onv.unsigned_representation);
    const char* end = begin + onv.K;
    while(begin != end) {
        os << std::bitset<CHAR_BIT>(*begin++);
    }
    return os;
}



/*
 *  SETTERS & GETTERS
 */

/**
 *  @set to a new representation
 */
void ONV::set_representation(size_t unsigned_representation) {
    this->unsigned_representation = unsigned_representation;
    update();
}

/**
 *  @return occupied orbital based on the electron index
 */
size_t ONV::get_occupied_orbital(size_t electron_index) {
    return occupation_indexes(electron_index);
}



// PUBLIC METHODS (PREVIOUS SPIN STRING FUNCTIONALITY)
/**
 *  @return if the index is occupied (i.e. 1) for the @param p-th spatial orbital, starting from 0
 *  @param p is the lexical index (i.e. read from right to left)
 */
bool ONV::isOccupied(size_t p) const{
    if (p > this->K-1) {
        throw std::invalid_argument("The index is out of the bitset bounds");
    }

    size_t operator_string = 1U << p;
    return this->unsigned_representation & operator_string;
}


/**
 *  @return if we can apply the annihilation operator (i.e. 1->0) for the @param p-th spatial orbital, starting from 0
 *  performs an annihilation @param p
 *
 *  !!! IMPORTANT: does not update the occupation array if required call "update()" !!!
 *  !!! IMPORTANT: does not update the occupation array if required call "update()" !!!
 */
bool ONV::annihilate(size_t p){

    if (this->isOccupied(p)) {
        size_t operator_string = 1U << p;
        this->unsigned_representation &= ~operator_string;
        return true;
    }else{
        return false;
    }
}

/**
 *  @return if we can apply the annihilation operator (i.e. 1->0) for the @param p-th spatial orbital, starting from 0
 *  @param p is the lexical index (i.e. read from right to left)
 *
 *  Furthermore, the @param sign is changed according to the sign change (+1 or -1) of the spin string after annihilation.
 *
 *  !!! IMPORTANT: does not update the occupation array if required call "update()" !!!
 *  !!! IMPORTANT: performs the annihilation in place !!!
 */
bool ONV::annihilate(size_t p, int& sign) {

    if (this->annihilate(p)) {  // we have to first check if we can annihilate before applying the phase factor
        sign *= this->operatorPhaseFactor(p);
        return true;
    } else {
        return false;
    }
}

/**
 * @return if we can apply the creation operator (i.e. 0->1) for the @param p-th spatial orbital, starting from 0
 *  performs a creation @param p
 *
 *  !!! IMPORTANT: does not update the occupation array if required call "update()" !!!
 *  !!! IMPORTANT: does not update the occupation array if required call "update()" !!!
 */
bool ONV::create(size_t p) {
    if (!this->isOccupied(p)) {
        size_t operator_string = 1U << p;
        this->unsigned_representation ^= operator_string;
        return true;
    } else {
        return false;
    }
}

/**
 *  @return if we can apply the creation operator (i.e. 0->1) for the @param p-th spatial orbital, starting from 0
 *  @param p is the lexical index (i.e. read from right to left)
 *
 *  Furthermore, the @param sign is changed according to the sign change (+1 or -1) of the spin string after annihilation.
 *
 *  !!! IMPORTANT: does not update the occupation array if required call "update()" !!!
 *  !!! IMPORTANT: performs the creation in place !!!
 */
bool ONV::create(size_t p, int& sign) {

    if (this->create(p)) {  // we have to first check if we can create before applying the phase factor
        sign *= this->operatorPhaseFactor(p);
        return true;
    } else {
        return false;
    }
}

/**
 *  @return the phase factor (+1 or -1) that arises by applying an annihilation or creation operator on orbital @param p, starting from 0 in reverse lexical ordering.
 *
 *  Let's say that there are m electrons in the orbitals up to p (not included). If m is even, the phase factor is (+1) and if m is odd, the phase factor is (-1), since electrons are fermions.
 */
int ONV::operatorPhaseFactor(size_t p) const {

    if (p == 0) {  // we can't give this to this->slice(0, 0)
        return 1;
    }
    size_t m = __builtin_popcountl(this->slice(0, p));  // count the number of set bits in the slice [0,p-1]

    if ( m % 2 == 0 ) {  // even number of electrons: phase factor (+1)
        return 1;
    } else {  // odd number of electrons: phase factor (-1)
        return -1;
    }
}


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
size_t ONV::slice(size_t index_start, size_t index_end) const {

    // First, do some checks
    if (index_end <= index_start) {
        throw std::invalid_argument("index_end should be larger than index_start.");
    }

    if (index_end > this->K + 1) {
        throw std::invalid_argument("The last slicing index index_end cannot be greater than the number of spatial orbitals K.");
    }

    // The union of these conditions also include the case that index_start > this->K


    // Shift bits to the right
    size_t u = this->unsigned_representation >> index_start;


    // Create the correct mask
    size_t mask_length = index_end - index_start;
    size_t mask = ((1U) << mask_length) - 1;


    // Use the mask
    return u & mask;
}



}  // namespace GQCG
