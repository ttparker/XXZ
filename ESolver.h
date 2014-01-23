class Sector
{
	public:
		Sector() {};

	private:
		std::vector<int> positions;			// which rows and columns
											// of matrix are in sector
		int multiplicity;					// size of this symmetry sector
        Eigen::MatrixXd sectorMat;          // sector operator
        static int fullMatrixSize;
        double lancTolerance;
        
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver; // DM eigensystem
		Eigen::VectorXd sectorEvals;		// eigenvalues
		int sectorColumnCounter;			// tracks which sector eigenvector
											// to fill into a matrix eigenvector

		Sector(const std::vector<int>& qNumList, int qNum,
			   const Eigen::MatrixXd& mat, double lancTolerance = 0.);
		Eigen::VectorXd fillOutEvec(bool takeLowest);

	friend class HamSolver;
	friend class DMSolver;
};

class HamSolver
{
	private:
		std::pair<Eigen::VectorXd, double> gState;

		HamSolver(const Eigen::MatrixXd& mat, const std::vector<int>& qNumList,
				  int targetQNum, double lancTolerance);
	
	friend class TheBlock;
	friend class EffectiveHamiltonian;
};

class DMSolver
{
	private:
		Eigen::MatrixXd highestEvecs;
		std::vector<int> highestEvecQNums;

		DMSolver(const Eigen::MatrixXd& mat, const std::vector<int>& qNumList,
				 int evecsToKeep);

	friend class TheBlock;
};
