#ifndef MAIN_H
#define MAIN_H
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <unsupported/Eigen/KroneckerProduct>
#include "d.h"

#define kp kroneckerProduct
#define Id(size) MatrixXd::Identity(size, size)
#define Id_d MatrixDd::Identity()   // one-site identity matrix

typedef Eigen::Matrix<double, d, d> MatrixDd;

#endif
