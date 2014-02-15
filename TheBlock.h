typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
    rmMatrixXd;

class EffectiveHamiltonian;

class TheBlock
{
	public:
        TheBlock(int m = 0,
				 const Eigen::MatrixXd& hS = Eigen::MatrixXd(),
				 const std::vector<Eigen::MatrixXd>& rhoBasisH2 
						= std::vector<Eigen::MatrixXd>(),
				 const std::vector<int>& qNumList = std::vector<int>());
		TheBlock(const Hamiltonian& ham, int mMaxIn);
        TheBlock nextBlock(const Hamiltonian& ham, bool exactDiag = true,
                           bool infiniteStage = true, int l = 0,
						   const TheBlock& compBlock = TheBlock(),
                           const TheBlock& beforeCompBlock = TheBlock());
                                                   // performs each DMRG step
        void reflectPredictedPsi();            // when you reach edge of system
        EffectiveHamiltonian createHSuperFinal(const Hamiltonian& ham,
                                               int skips) const;
					// HSuperFinal, mSFinal, qNumList, oneSiteQNums, targetQNum

	private:
		Eigen::MatrixXd hS;								// block Hamiltonian
		std::vector<Eigen::MatrixXd> rhoBasisH2;
									// density-matrix-basis coupling operators
		int m;								// number of states stored in block
		static rmMatrixXd psiGround;
        static bool firstfDMRGStep;         // slight abuse of nomenclature -
                                            // true during iDMRG as well
        std::vector<int> qNumList;			// tracks the conserved quantum
											// number of each row/column of hS
		static int mMax;				// max size of effective Hamiltonian
		Eigen::MatrixXd primeToRhoBasis;			// change-of-basis matrix

		Eigen::MatrixXd changeBasis(const Eigen::MatrixXd& mat) const;
				// represents operators in the basis of the new system block

	friend class EffectiveHamiltonian;
};
