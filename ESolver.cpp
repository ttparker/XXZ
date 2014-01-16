#include <set>
#include <map>
#include "d.h"
#include "main.h"
#include "Hamiltonian.h"
#include "TheBlock.h"
#include "ESolver.h"

using namespace Eigen;

Sector::Sector(const std::vector<int>& qNumList, int qNum, const MatrixXd& mat)
	: multiplicity(std::count(qNumList.begin(), qNumList.end(), qNum))
{
	positions.reserve(multiplicity);
	for(auto firstElement = qNumList.begin(), qNumListElement = firstElement,
		end = qNumList.end(); qNumListElement != end; qNumListElement++)
		if(*qNumListElement == qNum)
		   positions.push_back(qNumListElement - firstElement);
	MatrixXd sectorMat(multiplicity, multiplicity);			// sector operator
	int elmt = 0;
	for(int j : positions)			// fill sector matrix elements from matrix
		for(int i : positions)
			sectorMat(elmt++) = mat(i, j);
	SelfAdjointEigenSolver<MatrixXd> solver(sectorMat);
	sectorEvecs = solver.eigenvectors();
	sectorEvals = solver.eigenvalues();						// and eigenvalues
	sectorColumnCounter = multiplicity;		// start from the rightmost column
	fullMatrixSize = mat.rows();
};

VectorXd Sector::fillOutEvec(bool takeLowestIn)
{
	int sectorColumn = takeLowestIn ? 0 : --sectorColumnCounter;
	VectorXd evec = VectorXd::Zero(fullMatrixSize);
	for(int i = 0; i < multiplicity; i++)
		evec(positions[i]) = sectorEvecs(i, sectorColumn);
	return evec;
};

HamSolver::HamSolver(const Eigen::MatrixXd& mat, const std::vector<int>& qNumList,
					 int targetQNum)
{
	Sector::fullMatrixSize = mat.rows();
	Sector targetSector(qNumList, targetQNum, mat);
	lowestEval = targetSector.sectorEvals(0);
	psiGround = targetSector.fillOutEvec(true);
					// fill out full matrix eigenvector from stored sector one
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
		for(int i = 0, end = sectors[qNum].sectorEvals.size(); i < end; i++)
			indexedEvals.insert(		// add indexed eigenvalues to list
					std::pair<double, int>(sectors[qNum].sectorEvals(i), qNum));
	};
	highestEvecQNums.reserve(evecsToKeep);
	highestEvecs = MatrixXd::Zero(matSize, evecsToKeep);
    auto currentIndexedEval = indexedEvals.rbegin();
    for(int j = 0; j < evecsToKeep; j++)
	{
		int qNum = currentIndexedEval++ -> second;
		highestEvecQNums.push_back(qNum);
		highestEvecs.col(j) = sectors[qNum].fillOutEvec(false);
	};
};
