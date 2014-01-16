#ifndef INCLUDES_IN
#define INCLUDES_IN
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <unsupported/Eigen/KroneckerProduct>
#endif

#define kp kroneckerProduct
#define Id(size) MatrixXd::Identity(size, size)
#define Id_d MatrixDd::Identity()   // one-site identity matrix

typedef Eigen::Matrix<double, d, d> MatrixDd;
