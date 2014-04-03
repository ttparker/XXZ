#include "FreeFunctions.h"
#include "ESolver.h"

using namespace Eigen;

Hamiltonian TheBlock::ham;
rmMatrixXd TheBlock::psiGround;
int TheBlock::mMax;
bool TheBlock::firstfDMRGStep;

TheBlock::TheBlock(int m, const MatrixXd& hS,
                   const std::vector<MatrixXd>& rhoBasisH2,
                   const std::vector<int>& qNumList)
    : qNumList(qNumList), hS(hS), rhoBasisH2(rhoBasisH2), m(m) {};

TheBlock::TheBlock(const Hamiltonian& hamIn, int mMaxIn)
    : qNumList(ham.oneSiteQNums), hS(MatrixDd::Zero()), m(d)
{
    firstfDMRGStep = true;
    ham = hamIn;
    mMax = mMaxIn;
    rhoBasisH2.assign(ham.h2.begin(), ham.h2.begin() + indepCouplingOperators);
};

TheBlock TheBlock::nextBlock(const TheBlock& compBlock, int l, bool exactDiag,
                             bool infiniteStage,
                             const TheBlock& beforeCompBlock)
{
    std::vector<int> hSprimeQNumList	// add in quantum numbers of new site
        = vectorProductSum(qNumList, ham.oneSiteQNums);
    MatrixXd hSprime = kp(hS, Id_d)	+ ham.blockSiteJoin(rhoBasisH2);
                                                    // expanded system block
    std::vector<MatrixXd> tempRhoBasisH2;
    tempRhoBasisH2.reserve(indepCouplingOperators);
    int md = m * d;
    if(exactDiag)
    { // if near edge of system, no truncation necessary so skip DMRG algorithm
        for(auto op = ham.h2.begin(), end = op + indepCouplingOperators;
            op != end; op++)
            tempRhoBasisH2.push_back(kp(Id(m), *op));
        return TheBlock(md, hSprime, tempRhoBasisH2, hSprimeQNumList);
    };
    int compm = compBlock.m,
        compmd = compm * d;
    if(infiniteStage)
    {
        psiGround = VectorXd::Random(md * md);
        psiGround /= psiGround.norm();
    }
    else if(firstfDMRGStep)
    {
        psiGround = VectorXd::Random(md * compmd);
        psiGround /= psiGround.norm();
        firstfDMRGStep = false;
    }
    HamSolver hSuperSolver = (infiniteStage ?    // find superblock eigenstates
                              HamSolver(MatrixXd(kp(hSprime, Id(md))
                                                 + ham.siteSiteJoin(m, m)
                                                 + kp(Id(md), hSprime)),
                                        vectorProductSum(hSprimeQNumList,
                                                         hSprimeQNumList),
                                        ham.targetQNum * (l + 2) / ham.lSys * 2,
                                        psiGround) : // int automatically rounds down
                              HamSolver(MatrixXd(kp(hSprime, Id(compmd))
                                                 + ham.siteSiteJoin(m, compm)
                                                 + kp(Id(md), ham.blockSiteJoin(compBlock.rhoBasisH2)
                                                              + kp(compBlock.hS, Id_d))),
                                        vectorProductSum(hSprimeQNumList,
                                                         vectorProductSum(compBlock.qNumList,
                                                                          ham.oneSiteQNums)),
                                        ham.targetQNum, psiGround));
    psiGround = hSuperSolver.lowestEvec;                        // ground state
    psiGround.resize(md, infiniteStage ? md : compmd);
    DMSolver rhoSolver(psiGround * psiGround.adjoint(), hSprimeQNumList, mMax);
                                             // find density matrix eigenstates
    primeToRhoBasis = rhoSolver.highestEvecs; // construct change-of-basis matrix
    for(auto op = ham.h2.begin(), end = op + indepCouplingOperators; op != end;
        op++)
        tempRhoBasisH2.push_back(changeBasis(kp(Id(m), *op)));
    if(!infiniteStage)     // modify psiGround to predict the next ground state
    {
        for(int sPrimeIndex = 0; sPrimeIndex < md; sPrimeIndex++)
                    // transpose the environment block and right-hand free site
        {
            rmMatrixXd ePrime = psiGround.row(sPrimeIndex);
            ePrime.resize(compm, d);
            ePrime.transposeInPlace();
            ePrime.resize(1, compmd);
            psiGround.row(sPrimeIndex) = ePrime;
        };
        psiGround = primeToRhoBasis.adjoint() * psiGround; 
                                      // change the expanded system block basis
        psiGround.resize(mMax * d, compm);
        psiGround *= beforeCompBlock.primeToRhoBasis.transpose();
                                          // change the environment block basis
        psiGround.resize(mMax * d * beforeCompBlock.primeToRhoBasis.rows(), 1);
    };
    return TheBlock(mMax, changeBasis(hSprime), tempRhoBasisH2,
                    rhoSolver.highestEvecQNums);
                                // save expanded-block operators in new basis
};

void TheBlock::randomSeed(const TheBlock& compBlock)
{
    psiGround = VectorXd::Random(m * d * compBlock.m * d);
    psiGround /= psiGround.norm();
};

void TheBlock::reflectPredictedPsi()
{
    psiGround.resize(mMax * d, m * d);
    psiGround.transposeInPlace();
    psiGround.resize(mMax * d * m * d, 1);
};

EffectiveHamiltonian TheBlock::createHSuperFinal(const TheBlock& compBlock,
                                                 int skips) const
{
    int compm = compBlock.m;
    return EffectiveHamiltonian(qNumList, compBlock.qNumList, ham,
                                MatrixXd(kp(kp(hS, Id_d)
                                            + ham.blockSiteJoin(rhoBasisH2), Id(compm * d))
                                         + ham.siteSiteJoin(m, compm)
                                         + kp(Id(m * d), ham.blockSiteJoin(compBlock.rhoBasisH2)
                                                         + kp(compBlock.hS, Id_d))),
                                m, compm, skips);
};

MatrixXd TheBlock::changeBasis(const MatrixXd& mat) const
{
    return primeToRhoBasis.adjoint() * mat * primeToRhoBasis;
};
