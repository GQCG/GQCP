#ifndef GQCG_COMMON_HPP
#define GQCG_COMMON_HPP


#include <cstdlib>
#include <vector>

#include <Eigen/Dense>


namespace GQCG {


typedef std::vector<size_t> Vectoru;
typedef std::vector<Vectoru> Matrixu;
using VectorXs = Eigen::Matrix<size_t, Eigen::Dynamic, 1>;


}  // namespace GQCG


#endif  // GQCG_COMMON_HPP
