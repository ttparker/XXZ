#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include "main.h"

#define kp kroneckerProduct
#define Id(size) MatrixXd::Identity(size, size)

typedef std::vector<MatrixD_t, Eigen::aligned_allocator<MatrixD_t>> vecMatD_t;

class Hamiltonian
{
    public:
        Hamiltonian();
        void setParams(const std::vector<double>& couplingConstants,
                       int targetQNumIn, int lSysIn);
        
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    
    private:
        std::vector<int> oneSiteQNums;              // one-site quantum numbers
        std::vector<double> couplingConstants;
        int targetQNum,                              // targeted quantum number
            lSys;                                      // current system length
        vecMatD_t siteBasisH2;                 // site-basis coupling operators
        
        MatrixX_t blockSiteJoin(const std::vector<MatrixX_t>& rhoBasisH2) const,
                                           // appends free site to system block
                  siteSiteJoin(int m, int compm) const;
                                           // joins the two free sites together
    
    friend class TheBlock;
};

#endif
