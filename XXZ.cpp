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

void Hamiltonian::setParams(const std::vector<double>& couplingConstantsIn,
                            int targetQNumIn, int lSysIn)
{
    couplingConstants = couplingConstantsIn;
    targetQNum = targetQNumIn;
    lSys = lSysIn;
};

MatrixXd Hamiltonian::blockSiteJoin(const std::vector<MatrixXd>& rhoBasisH2) const
{
	return jz * kp(rhoBasisSigmaz, sigmaz)
		   + 2 * jxy * (kp(rhoBasisSigmaplus, sigmaminus)
						+ kp(rhoBasisSigmaplus.adjoint(), sigmaplus));
};

MatrixXd Hamiltonian::siteSiteJoin(int ml, int mlE) const
{
	return kp(Id(ml), jz * kp(kp(sigmaz, Id(mlE)), sigmaz)
                      + 2 * jxy * (kp(kp(sigmaplus, Id(mlE)), sigmaminus)
                                   + kp(kp(sigmaminus, Id(mlE)), sigmaplus)));
};
