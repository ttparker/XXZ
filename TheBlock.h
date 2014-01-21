typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
    rmMatrixXd;

class TheBlock
{
	public:
		TheBlock(int m = 0,
				 const Eigen::MatrixXd& hS = Eigen::MatrixXd(),
				 const std::vector<Eigen::MatrixXd>& rhoBasisH2 
						= std::vector<Eigen::MatrixXd>(),
				 const std::vector<int>& qNumList = std::vector<int>());
		TheBlock(const Hamiltonian& ham, int mMaxIn);
        TheBlock nextBlock(const Hamiltonian& ham, bool infiniteStage,
						   const TheBlock& compBlock, int l);	// performs each DMRG step
				// the third argument is the environment block in the fDMRG stage
		std::tuple<Eigen::MatrixXd, int, std::vector<int>, std::vector<int>, int>
			createHSuperFinal(const Hamiltonian& ham) const;
					// HSuperFinal, mSFinal, qNumList, oneSiteQNums, targetQNum

	private:
		Eigen::MatrixXd hS;								// block Hamiltonian
		std::vector<Eigen::MatrixXd> rhoBasisH2;
									// density-matrix-basis coupling operators
		int m;								// number of states stored in block
		std::vector<int> qNumList;			// tracks the conserved quantum
											// number of each row/column of hS
		static int mMax;				// max size of effective Hamiltonian
		Eigen::MatrixXd primeToRhoBasis;			// change-of-basis matrix

		Eigen::MatrixXd changeBasis(const Eigen::MatrixXd& mat) const;
				// represents operators in the basis of the new system block

	friend class EffectiveHamiltonian;
};
