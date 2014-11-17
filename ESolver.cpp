#include <map>
#include <set>
#include "GlobalPrecisionParameters.h"
#include "main.h"
#include "FreeFunctions.h"

using namespace Eigen;

std::vector<int> findTargetQNumPositions(const std::vector<int>& qNumList,
                                         int targetQNum)
{
    std::vector<int> positions; // which rows and columns of matrix are in sector
    for(auto firstElement = qNumList.begin(),
        qNumListElement = firstElement, end = qNumList.end();
        qNumListElement != end; qNumListElement++)
        if(*qNumListElement == targetQNum)
            positions.push_back(qNumListElement - firstElement);
    return positions;
};

MatrixX_t createCrossSectorBlock(const MatrixX_t& mat,
                                 const std::vector<int>& rowPositions,
                                 const std::vector<int>& colPositions)
                                         // pick out one block from full matrix
{
    MatrixX_t crossSectorBlock(rowPositions.size(), colPositions.size());
    int elmt = 0;
    for(int j : colPositions)
        for(int i : rowPositions)
            crossSectorBlock(elmt++) = mat(i, j);
    return crossSectorBlock;
};

VectorX_t filledOutEvec(const VectorX_t& sectorEvec, int fullMatrixSize,
                        const std::vector<int>& positions)
{
    VectorX_t longEvec = VectorX_t::Zero(fullMatrixSize);
    for(int i = 0, end = sectorEvec.size(); i < end; i++)
        longEvec(positions[i]) = sectorEvec(i);
    return longEvec;
};

HamSolver::HamSolver(const MatrixX_t& mat,
                     const std::vector<int>& hSprimeQNumList,
                     const std::vector<int>& hEprimeQNumList,
                     int targetQNum, rmMatrixX_t& bigSeed, double lancTolerance)
    : hSprimeQNumList(hSprimeQNumList), hEprimeQNumList(hEprimeQNumList),
      targetQNum(targetQNum)
{
    std::vector<int> positions
        = findTargetQNumPositions(vectorProductSum(hSprimeQNumList,
                                                   hEprimeQNumList),
                                  targetQNum);
    MatrixX_t sectorMat = createCrossSectorBlock(mat, positions, positions);
    int multiplicity = positions.size();
    rmMatrixX_t littleSeed(multiplicity, 1);
    for(int i = 0; i < multiplicity; i++)
        littleSeed(i) = bigSeed(positions[i]);
    littleSeed.normalize();
    lowestEval = lanczos(sectorMat, littleSeed, lancTolerance);
    lowestEvec = filledOutEvec(littleSeed, mat.rows(), positions);
};

DMSector::DMSector(const MatrixX_t& sectorMat, const std::vector<int>& positions)
    : sectorMat(sectorMat), positions(positions), multiplicity(positions.size()),
      sectorColumnCounter(multiplicity) {};

void DMSector::solveForAll()
{
    solver.compute(sectorMat);
};

VectorX_t DMSector::nextHighestEvec(int fullMatrixSize)
{
    return filledOutEvec(solver.eigenvectors().col(--sectorColumnCounter),
                         fullMatrixSize, positions);
};

DMSolver::DMSolver(const HamSolver hSuperSolver, int maxEvecsToKeep)
    : truncationError(0.)
{
    std::map<int, DMSector> sectors;        // key is the sector's quantum number
    std::multimap<double, int> indexedEvals;         // eigenvalue, then sector
    std::set<int> qNumSet(hSuperSolver.hSprimeQNumList.begin(),
                          hSuperSolver.hSprimeQNumList.end());
    for(int qNum : qNumSet)                 // make list of indexed eigenvalues
        if(std::count(hSuperSolver.hEprimeQNumList.begin(),
                      hSuperSolver.hEprimeQNumList.end(),
                      hSuperSolver.targetQNum - qNum))
        {
            std::vector<int> rowPositions
                = findTargetQNumPositions(hSuperSolver.hSprimeQNumList, qNum);
            MatrixX_t crossSectorBlock
                = createCrossSectorBlock(hSuperSolver.lowestEvec,
                                         rowPositions,
                                         findTargetQNumPositions(hSuperSolver
                                                                 .hEprimeQNumList,
                                                                 hSuperSolver.targetQNum
                                                                 - qNum));
            sectors.insert(sectors.end(),
                           std::make_pair(qNum,
                                          DMSector(crossSectorBlock
                                                   * crossSectorBlock.adjoint(),
                                                   rowPositions)));
                                                               // create sector
            sectors[qNum].solveForAll();
            for(int i = 0, end = sectors[qNum].multiplicity; i < end; i++)
                indexedEvals.insert(std::pair<double, int>
                                    (sectors[qNum].solver.eigenvalues()(i),
                                     qNum)); // add indexed eigenvalues to list
        };
    int matSize = hSuperSolver.lowestEvec.rows(),
        evecsToKeep;
    auto weight = indexedEvals.crbegin();
    if(matSize <= maxEvecsToKeep)
        evecsToKeep = matSize;
    else
    {
        evecsToKeep = maxEvecsToKeep;
        std::advance(weight, maxEvecsToKeep - 1);
        for(; evecsToKeep >= 1
              && (weight -> first == 0
                  || (weight -> first - std::next(weight, 1) -> first)
                      / std::abs(weight -> first) < degenerateDMCutoff);
            evecsToKeep--, weight--);
              // find the the max number of eigenvectors to keep that do not
              // terminate inside a degenerate eigenspace of the density matrix
        if(evecsToKeep == 0)
        {
            std::cerr << "More than mMax highest-weighted density-matrix "
                      << "eigenvectors are degenerate." << std::endl;
            exit(EXIT_FAILURE);
        }
        else if(evecsToKeep != maxEvecsToKeep)
            std::cout << "Warning: mMax truncation ends in a degenerate DM "
                      << "eigenspace, lowering cutoff to " << evecsToKeep
                      << " states." << std::endl;
    };
    weight++;                    // now points to first truncated DM eigenvalue
    for(auto end = indexedEvals.crend(); weight != end; weight++)
        truncationError += weight -> first;
    highestEvecQNums.reserve(evecsToKeep);
    highestEvecs = MatrixX_t::Zero(matSize, evecsToKeep);
    auto currentIndexedEval = indexedEvals.crbegin();
    for(int j = 0; j < evecsToKeep; j++)
    {
        int qNum = currentIndexedEval++ -> second;
        highestEvecQNums.push_back(qNum);
        highestEvecs.col(j) = sectors[qNum].nextHighestEvec(matSize);
    };
};
