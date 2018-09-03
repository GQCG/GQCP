#ifndef GQCG_UNITS_HPP
#define GQCG_UNITS_HPP


namespace GQCG {
namespace units {


namespace constants {

struct CODATA2014 {
    constexpr static double angstrom_per_bohr = 0.52917721067;
    constexpr static double bohr_per_angstrom = 1 / angstrom_per_bohr;
};

}  // namespace constants



/**
 *  Given a @param value_in_bohr, @return the value in Angstrom
 */
inline double bohr_to_angstrom(double value_in_bohr) { return value_in_bohr * GQCG::units::constants::CODATA2014::angstrom_per_bohr; }

/**
 *  Given a @param value_in_angstrom, @return the value in Bohr (a.u.)
 */
inline double angstrom_to_bohr(double value_in_angstrom) { return value_in_angstrom * GQCG::units::constants::CODATA2014::bohr_per_angstrom; }



}  // namespace units
}  // namespace GQCG


#endif  // GQCG_UNITS_HPP
