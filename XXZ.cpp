#include "Hamiltonian.h"

#define jxy couplingConstants[0]
#define jz couplingConstants[1]
#define sigmaplus siteBasisH2[0]
#define sigmaz siteBasisH2[1]
#define sigmaminus siteBasisH2[2]
#define rhoBasisSigmaplus rhoBasisH2[0]
#define rhoBasisSigmaz rhoBasisH2[1]

using namespace Eigen;

Hamiltonian::Hamiltonian() : oneSiteQNums({1, -1})
{
    siteBasisH2.resize(3);
    sigmaplus << 0., 1.,
                 0., 0.;
    sigmaminus << 0., 0.,
                  1., 0.;
    sigmaz << 1.,  0.,
              0., -1.;                                 // define Pauli matrices
};

void Hamiltonian::setParams(const std::vector<double>& couplingConstantsIn,
                            int targetQNumIn, int lSysIn)
{
    couplingConstants = couplingConstantsIn;
    targetQNum = targetQNumIn;
    lSys = lSysIn;
};

MatrixX_t Hamiltonian::blockSiteJoin(const std::vector<MatrixX_t>& rhoBasisH2)
    const
{
    MatrixX_t plusMinus = kp(rhoBasisSigmaplus, sigmaminus);
    return jz * kp(rhoBasisSigmaz, sigmaz) + 2 * jxy * (plusMinus
                                                        + plusMinus.adjoint());
};

MatrixX_t Hamiltonian::siteSiteJoin(int m, int compm) const
{
    MatrixX_t plusMinus = kp(kp(sigmaplus, Id(compm)), sigmaminus);
    return kp(Id(m), jz * kp(kp(sigmaz, Id(compm)), sigmaz)
                     + 2 * jxy * (plusMinus + plusMinus.adjoint()));
};
