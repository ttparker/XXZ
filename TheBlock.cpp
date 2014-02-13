#include <tuple>
#include "d.h"
#include "main.h"
#include "Hamiltonian.h"
#include "TheBlock.h"
#include "EffectiveHamiltonian.h"
#include "FreeFunctions.h"
#include "ESolver.h"

using namespace Eigen;

int TheBlock::mMax;
double TheBlock::lancTolerance;

TheBlock::TheBlock(int m, const MatrixXd& hS,
				   const std::vector<MatrixXd>& rhoBasisH2,
				   const std::vector<int>& qNumList)
	: hS(hS), rhoBasisH2(rhoBasisH2), m(m), qNumList(qNumList) {};

TheBlock::TheBlock(const Hamiltonian& ham, int mMaxIn)
	: hS(MatrixDd::Zero()), m(d), qNumList(ham.oneSiteQNums)
{
	mMax = mMaxIn;
	rhoBasisH2.assign(ham.h2.begin(),
					  ham.h2.begin() + ham.couplingConstants.size());
};

TheBlock TheBlock::nextBlock(const Hamiltonian& ham, bool exactDiag,
                             bool infiniteStage, int l,
                             const TheBlock& compBlock)
{
	std::vector<int> hSprimeQNumList	// add in quantum numbers of new site
		= vectorProductSum(qNumList, ham.oneSiteQNums);
    MatrixXd hSprime = kp(hS, Id_d)	+ ham.blockSiteJoin(rhoBasisH2);
													// expanded system block
	std::vector<MatrixXd> tempRhoBasisH2;
	int indepCouplingOperators = ham.couplingConstants.size();
	tempRhoBasisH2.reserve(indepCouplingOperators);
	int md = m * d;
	if(exactDiag)
	{ // if near edge of system, no truncation necessary so skip DMRG algorithm
        for(auto op = ham.h2.begin(), end = ham.h2.begin() + indepCouplingOperators;
			op != end; op++)
		    tempRhoBasisH2.push_back(kp(Id(m), *op));
		return TheBlock(md, hSprime, tempRhoBasisH2, hSprimeQNumList);
	};
    int compmd = compBlock.m * d;
	HamSolver hSuperSolver = (infiniteStage ?		// find superblock eigenstates
							  HamSolver(MatrixXd(kp(hSprime, Id(md))
												 + ham.siteSiteJoin(m, m)
												 + kp(Id(md), hSprime)),
										vectorProductSum(hSprimeQNumList,
														hSprimeQNumList),
										ham.targetQNum * (l + 2) / ham.lSys * 2,
                                        lancTolerance) :
											// int automatically rounds down
							  HamSolver(MatrixXd(kp(hSprime, Id(compmd))
												 + ham.siteSiteJoin(m, compBlock.m)
												 + kp(Id(md),
													  ham.blockSiteJoin(compBlock.rhoBasisH2))
												 + kp(kp(Id(md), compBlock.hS),
													  Id_d)),
										vectorProductSum(hSprimeQNumList,
														 vectorProductSum(compBlock.qNumList,
																		  ham.oneSiteQNums)),
										ham.targetQNum,
                                        lancTolerance));
	rmMatrixXd psiGround = hSuperSolver.gState.first;				// ground state
    psiGround.resize(md, infiniteStage ? md : compmd);
	DMSolver rhoSolver(psiGround * psiGround.adjoint(), hSprimeQNumList, mMax);
											// find density matrix eigenstates
	primeToRhoBasis = rhoSolver.highestEvecs;	// construct change-of-basis matrix
	for(auto op = ham.h2.begin(), end = ham.h2.begin() + indepCouplingOperators;
        op != end; op++)
		tempRhoBasisH2.push_back(changeBasis(kp(Id(m), *op)));
	return TheBlock(mMax, changeBasis(hSprime), tempRhoBasisH2,
					rhoSolver.highestEvecQNums);
								// save expanded-block operators in new basis
};

EffectiveHamiltonian TheBlock::createHSuperFinal(const Hamiltonian& ham,
                                                 double lancTolerance,
                                                 int skips) const
{
	return EffectiveHamiltonian(qNumList, ham,
                                MatrixXd(kp(hS, Id(d * m * d))
								+ kp(ham.blockSiteJoin(rhoBasisH2), Id(m * d))
								+ ham.siteSiteJoin(m, m)
								+ kp(Id(m * d), ham.blockSiteJoin(rhoBasisH2))
								+ kp(kp(Id(m * d), hS), Id_d)),
                                lancTolerance, m, skips);
};

MatrixXd TheBlock::changeBasis(const MatrixXd& mat) const
{
	return primeToRhoBasis.adjoint() * mat * primeToRhoBasis;
};
