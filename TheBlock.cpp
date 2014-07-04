#include "FreeFunctions.h"
#include "ESolver.h"

using namespace Eigen;

TheBlock::TheBlock(int m, const std::vector<int>& qNumList, const MatrixX_t& hS,
                   const std::vector<MatrixX_t>& rhoBasisH2, int l)
    : m(m), qNumList(qNumList), hS(hS), rhoBasisH2(rhoBasisH2), l(l) {};

TheBlock::TheBlock(const Hamiltonian& ham)
    : m(d), qNumList(ham.oneSiteQNums), hS(MatrixD_t::Zero()), l(0)
{
    rhoBasisH2.assign(ham.h2.begin(), ham.h2.begin() + indepCouplingOperators);
};

TheBlock TheBlock::nextBlock(const stepData& data, rmMatrixX_t& psiGround)
{
    MatrixX_t hSprime = kp(hS, Id_d)
                        + data.ham.blockSiteJoin(rhoBasisH2);
                                                       // expanded system block
    std::vector<int> hSprimeQNumList = vectorProductSum(qNumList,
                                                        data.ham.oneSiteQNums);
                                          // add in quantum numbers of new site
    std::vector<MatrixX_t> tempRhoBasisH2;
    tempRhoBasisH2.reserve(indepCouplingOperators);
    int md = m * d;
    if(data.exactDiag)
      // if near edge of system, no truncation necessary so skip DMRG algorithm
    {
        for(auto op = data.ham.h2.begin(), end = op + indepCouplingOperators;
            op != end; op++)
            tempRhoBasisH2.push_back(kp(Id(m), *op));
        return TheBlock(md, hSprimeQNumList, hSprime, tempRhoBasisH2, l + 1);
    };
    int compm = data.compBlock -> m,
        compmd = compm * d;
    MatrixX_t hEprime = (data.infiniteStage ?
                         hSprime :
                         kp(data.compBlock -> hS, Id_d)
                         + data.ham.blockSiteJoin(data.compBlock -> rhoBasisH2));
                                                  // expanded environment block
    std::vector<int> hEprimeQNumList = (data.infiniteStage ?
                                        hSprimeQNumList :
                                        vectorProductSum(data.compBlock
                                                         -> qNumList,
                                                         data.ham.oneSiteQNums));
    int scaledTargetQNum = (data.infiniteStage ?
                            data.ham.targetQNum * (l + 2) / data.ham.lSys * 2 :
                                               // int automatically rounds down
                            data.ham.targetQNum);
                         // during iDMRG stage, targets correct quantum number
                         // per unit site by scaling to fit current system size
                         // - note: this will change if d != 2
    HamSolver hSuperSolver(kp(hSprime, Id(compmd))
                           + data.ham.siteSiteJoin(m, compm)
                           + kp(Id(md), hEprime),
                           vectorProductSum(hSprimeQNumList, hEprimeQNumList),
                           scaledTargetQNum, psiGround, data.lancTolerance);
                                                 // find superblock eigenstates
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
            rmMatrixX_t ePrime = psiGround.row(sPrimeIndex);
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
                                            const rmMatrixX_t& psiGround,
                                            int skips) const
{
    MatrixX_t hSprime = kp(hS, Id_d)
                        + data.ham.blockSiteJoin(rhoBasisH2);
                                                       // expanded system block
    int compm = data.compBlock -> m;
    MatrixX_t hEprime = kp(data.compBlock -> hS, Id_d)
                        + data.ham.blockSiteJoin(data.compBlock -> rhoBasisH2);
                                                  // expanded environment block
    return FinalSuperblock(kp(hSprime, Id(compm * d))
                           + data.ham.siteSiteJoin(m, compm)
                           + kp(Id(m * d), hEprime),
                           qNumList, data.compBlock -> qNumList, data,
                           psiGround, m, compm, skips);
};

obsMatrixX_t TheBlock::obsChangeBasis(const obsMatrixX_t& mat) const
{
    return primeToRhoBasis.adjoint() * mat * primeToRhoBasis;
};

MatrixX_t TheBlock::changeBasis(const MatrixX_t& mat) const
{
    return primeToRhoBasis.adjoint() * mat * primeToRhoBasis;
};
