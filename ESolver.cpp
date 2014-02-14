#include <map>
#include <set>
#include "d.h"
#include "main.h"
#include "Hamiltonian.h"
#include "TheBlock.h"
#include "ESolver.h"

using namespace Eigen;

double Sector::lancTolerance;
int Sector::fullMatrixSize;

Sector::Sector(const std::vector<int>& qNumList, int qNum, const MatrixXd& mat)
	: multiplicity(std::count(qNumList.begin(), qNumList.end(), qNum)),
	  sectorMat(MatrixXd(multiplicity, multiplicity)),
	  sectorColumnCounter(multiplicity)
{
	positions.reserve(multiplicity);
	for(auto firstElement = qNumList.begin(), qNumListElement = firstElement,
		end = qNumList.end(); qNumListElement != end; qNumListElement++)
		if(*qNumListElement == qNum)
		   positions.push_back(qNumListElement - firstElement);
	int elmt = 0;
	for(int j : positions)			// fill sector matrix elements from matrix
		for(int i : positions)
			sectorMat(elmt++) = mat(i, j);
};

VectorXd Sector::filledOutEvec(VectorXd sectorEvec) const
{
	VectorXd longEvec = VectorXd::Zero(fullMatrixSize);
	for(int i = 0; i < multiplicity; i++)
		longEvec(positions[i]) = sectorEvec(i);
	return longEvec;
};

double Sector::solveForLowest(VectorXd& lowestEvec)
{
    VectorXd seed = VectorXd::Random(multiplicity);
    seed /= seed.norm();
    double lowestEval = lanczos(sectorMat, seed);
    lowestEvec = filledOutEvec(seed);
    return lowestEval;
};

void Sector::solveForAll()
{
    solver.compute(sectorMat);
};

Eigen::VectorXd Sector::nextHighestEvec()
{
    return filledOutEvec(solver.eigenvectors().col(--sectorColumnCounter));
};

HamSolver::HamSolver(const Eigen::MatrixXd& mat, const std::vector<int>& qNumList,
					 int targetQNum)
{
	Sector::fullMatrixSize = mat.rows();
	Sector targetSector(qNumList, targetQNum, mat);
    lowestEval = targetSector.solveForLowest(lowestEvec);
};

DMSolver::DMSolver(const Eigen::MatrixXd& mat, const std::vector<int>& qNumList,
				   int evecsToKeep)
{
	int matSize = mat.rows();
	Sector::fullMatrixSize = matSize;
	std::map<int, Sector> sectors;			// key is the quantum number
	std::map<double, int> indexedEvals;		// eigenvalue, then sector
	std::set<int> qNumSet(qNumList.begin(), qNumList.end());
	for(int qNum : qNumSet)				// make list of indexed eigenvalues
	{
		sectors.insert(sectors.end(),					// create sector
					   std::make_pair(qNum, Sector(qNumList, qNum, mat)));
        sectors[qNum].solveForAll();
		for(int i = 0, end = sectors[qNum].multiplicity; i < end; i++)
			indexedEvals.insert(std::pair<double, int>
                                (sectors[qNum].solver.eigenvalues()(i), qNum));
                                            // add indexed eigenvalues to list
	};
	highestEvecQNums.reserve(evecsToKeep);
	highestEvecs = MatrixXd::Zero(matSize, evecsToKeep);
    auto currentIndexedEval = indexedEvals.rbegin();
    for(int j = 0; j < evecsToKeep; j++)
	{
		int qNum = currentIndexedEval++ -> second;
		highestEvecQNums.push_back(qNum);
		highestEvecs.col(j) = sectors[qNum].nextHighestEvec();
	};
};
