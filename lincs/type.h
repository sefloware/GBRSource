#ifndef TYPE_H
#define TYPE_H

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/lagged_fibonacci.hpp>
#include <Eigen/Eigen>
#include <Eigen/Sparse>

using namespace Eigen;
typedef Matrix<Vector3d,Dynamic,1> VectorXV3d;
typedef Matrix< bool, 3, 1 > Vector3b;
//typedef boost::mt19937 base_generator_type;
typedef boost::lagged_fibonacci44497 base_generator_type;

#define M_PI    3.14159265358979323846 //pi
#define M_KB    1.3806505E-23   //The boltzmann constant

#endif // TYPE_H
