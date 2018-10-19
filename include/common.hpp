#ifndef GQCP_COMMON_HPP
#define GQCP_COMMON_HPP


#include <cstdlib>
#include <vector>

#include <Eigen/Dense>


namespace GQCP {


typedef std::vector<size_t> Vectoru;
typedef std::vector<Vectoru> Matrixu;
using VectorXs = Eigen::Matrix<size_t, Eigen::Dynamic, 1>;


}  // namespace GQCP


#endif  // GQCP_COMMON_HPP
