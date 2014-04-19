#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include "main.h"

#define kp kroneckerProduct
#define Id(size) MatrixXd::Identity(size, size)

class Hamiltonian
{
    public:
        std::vector<int> oneSiteQNums;              // one-site quantum numbers
        int targetQNum,                              // targeted quantum number
            lSys;                                      // current system length
        
        Hamiltonian();
        void setParams(const std::vector<double>& couplingConstants,
                       int targetQNumIn, int lSysIn);
        
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    
    private:
        std::vector<double> couplingConstants;
        std::vector<MatrixDd, Eigen::aligned_allocator<MatrixDd>> h2;
                                               // site-basis coupling operators
        
        Eigen::MatrixXd
            blockSiteJoin(const std::vector<Eigen::MatrixXd>& rhoBasisH2) const,
                                           // appends free site to system block
            siteSiteJoin(int m1, int m2) const;
                                           // joins the two free sites together
    
    friend class TheBlock;
};

#endif
