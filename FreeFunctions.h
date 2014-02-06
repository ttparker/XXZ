std::vector<int> vectorProductSum(const std::vector<int>& first,
								  const std::vector<int>& second);
			// takes the tensor product of two blocks' lists of quantum numbers
void oneSiteExpValues(const MatrixDd& oneSiteOp, int rangeOfObservables,
					  int currentLSys, EffectiveHamiltonian& hSuperFinal,
					  std::vector<TheBlock>& blocks, std::ofstream& fileout);
void twoSiteExpValues(const MatrixDd& firstTwoSiteOp,
					  const MatrixDd& secondTwoSiteOp, int rangeOfObservables,
					  int currentLSys, EffectiveHamiltonian& hSuperFinal,
					  std::vector<TheBlock>& blocks, std::ofstream& fileout);
void modifyHamParams(int trial = 0);
