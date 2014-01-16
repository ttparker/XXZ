class Hamiltonian
{
	public:
		int lSys;									// current system length

		Hamiltonian(int lSys, const std::vector<double>& couplingConstants,
					int targetQNum);
 
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  
	private:
		std::vector<double> couplingConstants;
		std::vector<MatrixDd, Eigen::aligned_allocator<MatrixDd>> h2;
                                            // site-basis coupling operators
		int targetQNum;				// targeted average magnetization per site
		std::vector<int> oneSiteQNums;

		Eigen::MatrixXd
			blockSiteJoin(const std::vector<Eigen::MatrixXd>& rhoBasisH2) const,
										// appends free site to system block
			siteSiteJoin(int m1, int m2) const;
										// joins the two free sites together

	friend class TheBlock;
    friend void modifyHamParams(int trial);
};
