#include "FreeFunctions.h"
#include "ESolver.h"

using namespace Eigen;

TheBlock::TheBlock(int m, const std::vector<int>& qNumList, const MatrixXd& hS,
                   const std::vector<MatrixXd>& rhoBasisH2, int l)
    : m(m), qNumList(qNumList), hS(hS), rhoBasisH2(rhoBasisH2), l(l) {};

TheBlock::TheBlock(const Hamiltonian& ham)
    : m(d), qNumList(ham.oneSiteQNums), hS(MatrixDd::Zero()), l(0)
{
    rhoBasisH2.assign(ham.h2.begin(), ham.h2.begin() + indepCouplingOperators);
};

TheBlock TheBlock::nextBlock(const stepData& data, rmMatrixXd& psiGround)
{
    std::vector<int> hSprimeQNumList	  // add in quantum numbers of new site
        = vectorProductSum(qNumList, data.ham.oneSiteQNums);
    MatrixXd hSprime = kp(hS, Id_d)	+ data.ham.blockSiteJoin(rhoBasisH2);
                                                       // expanded system block
    std::vector<MatrixXd> tempRhoBasisH2;
    tempRhoBasisH2.reserve(indepCouplingOperators);
    int md = m * d;
    if(data.exactDiag)
    { // if near edge of system, no truncation necessary so skip DMRG algorithm
        for(auto op = data.ham.h2.begin(), end = op + indepCouplingOperators;
            op != end; op++)
            tempRhoBasisH2.push_back(kp(Id(m), *op));
        return TheBlock(md, hSprimeQNumList, hSprime, tempRhoBasisH2, l + 1);
    };
    int compm = data.compBlock -> m,
        compmd = compm * d;
    HamSolver hSuperSolver = (data.infiniteStage ? // find superblock eigenstates
                              HamSolver(MatrixXd(kp(hSprime, Id(md))
                                                 + data.ham.siteSiteJoin(m, m)
                                                 + kp(Id(md), hSprime)),
                                        vectorProductSum(hSprimeQNumList,
                                                         hSprimeQNumList),
                                        data.ham.targetQNum * (l + 2) / data.ham.lSys * 2,
                                        psiGround, data.lancTolerance) :
                                               // int automatically rounds down
                              HamSolver(MatrixXd(kp(hSprime, Id(compmd))
                                                 + data.ham.siteSiteJoin(m, compm)
                                                 + kp(Id(md), data.ham.blockSiteJoin(data.compBlock -> rhoBasisH2)
                                                              + kp(data.compBlock -> hS, Id_d))),
                                        vectorProductSum(hSprimeQNumList,
                                                         vectorProductSum(data.compBlock -> qNumList,
                                                                          data.ham.oneSiteQNums)),
                                        data.ham.targetQNum, psiGround,
                                        data.lancTolerance));
    psiGround = hSuperSolver.lowestEvec;                        // ground state
    psiGround.resize(md, compmd);
    DMSolver rhoSolver(psiGround * psiGround.adjoint(), hSprimeQNumList, data.mMax);
                                             // find density matrix eigenstates
    primeToRhoBasis = rhoSolver.highestEvecs; // construct change-of-basis matrix
    for(auto op = data.ham.h2.begin(), end = op + indepCouplingOperators;
        op != end; op++)
        tempRhoBasisH2.push_back(changeBasis(kp(Id(m), *op)));
    if(!data.infiniteStage)     // modify psiGround to predict the next ground state
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
        psiGround.resize(data.mMax * d, compm);
        psiGround *= data.beforeCompBlock -> primeToRhoBasis.transpose();
                                          // change the environment block basis
        psiGround.resize(data.mMax * d
                         * data.beforeCompBlock -> primeToRhoBasis.rows(), 1);
    };
    return TheBlock(data.mMax, rhoSolver.highestEvecQNums, changeBasis(hSprime),
                    tempRhoBasisH2, l + 1);
                                  // save expanded-block operators in new basis
};

FinalSuperblock TheBlock::createHSuperFinal(const stepData& data,
                                            const rmMatrixXd& psiGround,
                                            int skips) const
{
    int compm = data.compBlock -> m;
    return FinalSuperblock(MatrixXd(kp(kp(hS, Id_d)
                                       + data.ham.blockSiteJoin(rhoBasisH2),
                                       Id(compm * d))
                                    + data.ham.siteSiteJoin(m, compm)
                                    + kp(Id(m * d), data.ham.blockSiteJoin(data.compBlock -> rhoBasisH2)
                                                    + kp(data.compBlock -> hS, Id_d))),
                           qNumList, data.compBlock -> qNumList, data,
                           psiGround, m, compm, skips);
};

MatrixXd TheBlock::changeBasis(const MatrixXd& mat) const
{
    return primeToRhoBasis.adjoint() * mat * primeToRhoBasis;
};
