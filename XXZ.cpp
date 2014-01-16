#include "d.h"
#include "main.h"
#include "Hamiltonian.h"

#define jxy couplingConstants[0]
#define jz couplingConstants[1]
#define sigmaplus h2[0]
#define sigmaz h2[1]
#define sigmaminus h2[2]
#define rhoBasisSigmaplus rhoBasisH2[0]
#define rhoBasisSigmaz rhoBasisH2[1]

using namespace Eigen;

Hamiltonian::Hamiltonian(int lSys, const std::vector<double>& couplingConstants,
						 int targetQNum)
	: lSys(lSys), couplingConstants(couplingConstants), targetQNum(targetQNum)
{
	h2.resize(3);
	sigmaplus << 0., 1.,
				 0., 0.;
	sigmaminus << 0., 0.,
				  1., 0.;
	sigmaz << 1., 0.,
			  0., -1.;								 // define Pauli matrices
	oneSiteQNums.reserve(2);
	oneSiteQNums.push_back(1);
	oneSiteQNums.push_back(-1);
};

void Hamiltonian::modifyParams(int trial)
{
	trial;
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
