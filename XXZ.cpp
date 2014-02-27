#include "Hamiltonian.h"

#define jxy couplingConstants[0]
#define jz couplingConstants[1]
#define sigmaplus h2[0]
#define sigmaz h2[1]
#define sigmaminus h2[2]
#define rhoBasisSigmaplus rhoBasisH2[0]
#define rhoBasisSigmaz rhoBasisH2[1]

using namespace Eigen;

Hamiltonian::Hamiltonian() : oneSiteQNums({1, -1})
{
	h2.resize(3);
	sigmaplus << 0., 1.,
				 0., 0.;
	sigmaminus << 0., 0.,
				  1., 0.;
	sigmaz << 1., 0.,
			  0., -1.;								 // define Pauli matrices
};

void Hamiltonian::setParams(int lSysIn,
                            const std::vector<double>& couplingConstantsIn,
                            int targetQNumIn)
{
    lSys = lSysIn;
    couplingConstants = couplingConstantsIn;
    targetQNum = targetQNumIn;
};

MatrixXd Hamiltonian::blockSiteJoin(const std::vector<MatrixXd>& rhoBasisH2) const
{
	return jz * kp(rhoBasisSigmaz, sigmaz)
		   + 2 * jxy * (kp(rhoBasisSigmaplus, sigmaminus)
						+ kp(rhoBasisSigmaplus.adjoint(), sigmaplus));
};

MatrixXd Hamiltonian::siteSiteJoin(int m1, int m2) const
{
	return jz * kp(kp(Id(m1), sigmaz), kp(Id(m2), sigmaz))
		   + 2 * jxy * (kp(kp(Id(m1), sigmaplus), kp(Id(m2), sigmaminus))
						+ kp(kp(Id(m1), sigmaminus), kp(Id(m2), sigmaplus)));
};
