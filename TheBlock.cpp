#include "FreeFunctions.h"
#include "ESolver.h"

using namespace Eigen;

rmMatrixXd TheBlock::psiGround;
bool TheBlock::firstfDMRGStep;
int TheBlock::mMax;

TheBlock::TheBlock(int m, const MatrixXd& hS,
				   const std::vector<MatrixXd>& rhoBasisH2,
				   const std::vector<int>& qNumList)
	: hS(hS), rhoBasisH2(rhoBasisH2), m(m), qNumList(qNumList) {};

TheBlock::TheBlock(const Hamiltonian& ham, int mMaxIn)
	: hS(MatrixDd::Zero()), m(d), qNumList(ham.oneSiteQNums)
{
    firstfDMRGStep = true;
	mMax = mMaxIn;
	rhoBasisH2.assign(ham.h2.begin(),
					  ham.h2.begin() + ham.couplingConstants.size());
};

TheBlock TheBlock::nextBlock(const Hamiltonian& ham, bool exactDiag,
                             bool infiniteStage, int l,
                             const TheBlock& compBlock,
                             const TheBlock& beforeCompBlock)
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
    VectorXd seed;
    if(infiniteStage)
    {
        seed = VectorXd::Random(md * md);
        seed /= seed.norm();
    }
    else if(firstfDMRGStep)
    {
        seed = VectorXd::Random(md * compmd);
        seed /= seed.norm();
        firstfDMRGStep = false;
    }
    else
        seed = psiGround;
	HamSolver hSuperSolver = (infiniteStage ?		// find superblock eigenstates
							  HamSolver(MatrixXd(kp(hSprime, Id(md))
												 + ham.siteSiteJoin(m, m)
												 + kp(Id(md), hSprime)),
										vectorProductSum(hSprimeQNumList,
														hSprimeQNumList),
										ham.targetQNum * (l + 2) / ham.lSys * 2,
                                        seed) :	// int automatically rounds down
							  HamSolver(MatrixXd(kp(hSprime, Id(compmd))
												 + ham.siteSiteJoin(m, compBlock.m)
												 + kp(Id(md),
													  ham.blockSiteJoin(compBlock.rhoBasisH2))
												 + kp(kp(Id(md), compBlock.hS),
													  Id_d)),
										vectorProductSum(hSprimeQNumList,
														 vectorProductSum(compBlock.qNumList,
																		  ham.oneSiteQNums)),
										ham.targetQNum, seed));
	psiGround = hSuperSolver.lowestEvec;				// ground state
    psiGround.resize(md, infiniteStage ? md : compmd);
	DMSolver rhoSolver(psiGround * psiGround.adjoint(), hSprimeQNumList, mMax);
											// find density matrix eigenstates
	primeToRhoBasis = rhoSolver.highestEvecs;	// construct change-of-basis matrix
	for(auto op = ham.h2.begin(), end = ham.h2.begin() + indepCouplingOperators;
        op != end; op++)
		tempRhoBasisH2.push_back(changeBasis(kp(Id(m), *op)));
    if(!infiniteStage)     // modify psiGround to predict the next ground state
    {
        for(int sPrimeIndex = 0; sPrimeIndex < md; sPrimeIndex++)
                    // transpose the environment block and right-hand free site
        {
            rmMatrixXd ePrime = psiGround.row(sPrimeIndex);
            ePrime.resize(compBlock.m, d);
            ePrime.transposeInPlace();
            ePrime.resize(1, compmd);
            psiGround.row(sPrimeIndex) = ePrime;
        };
        psiGround = primeToRhoBasis.adjoint() * psiGround; 
                                      // change the expanded system block basis
        psiGround.resize(mMax * d, compBlock.m);
        psiGround *= beforeCompBlock.primeToRhoBasis.transpose();
                                          // change the environment block basis
        psiGround.resize(mMax * d * beforeCompBlock.primeToRhoBasis.rows(), 1);
    };
	return TheBlock(mMax, changeBasis(hSprime), tempRhoBasisH2,
					rhoSolver.highestEvecQNums);
								// save expanded-block operators in new basis
};

void TheBlock::reflectPredictedPsi()
{
    psiGround.resize(mMax * d, m * d);
    psiGround.transposeInPlace();
    psiGround.resize(mMax * d * m * d, 1);
};

EffectiveHamiltonian TheBlock::createHSuperFinal(const Hamiltonian& ham,
                                                 const TheBlock& compBlock,
                                                 int skips) const
{
	return EffectiveHamiltonian(qNumList, ham,
                                MatrixXd(kp(hS, Id(d * compBlock.m * d))
								+ kp(ham.blockSiteJoin(rhoBasisH2), Id(compBlock.m * d))
								+ ham.siteSiteJoin(m, compBlock.m)
								+ kp(Id(m * d), ham.blockSiteJoin(compBlock.rhoBasisH2))
								+ kp(kp(Id(m * d), compBlock.hS), Id_d)),
                                m, skips);
};

MatrixXd TheBlock::changeBasis(const MatrixXd& mat) const
{
	return primeToRhoBasis.adjoint() * mat * primeToRhoBasis;
};
