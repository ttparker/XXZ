#ifndef FREEFUNCTIONS_H
#define FREEFUNCTIONS_H

#include "EffectiveHamiltonian.h"

std::vector<int> vectorProductSum(const std::vector<int>& first,
                                  const std::vector<int>& second);
            // takes the tensor product of two blocks' lists of quantum numbers
void oneSiteExpValues(const MatrixDd& oneSiteOp, int rangeOfObservables,
                      int currentLSys, EffectiveHamiltonian& hSuperFinal,
                      std::vector<TheBlock>& leftBlocks,
                      std::vector<TheBlock>& rightBlocks,
                      std::ofstream& fileout);
void twoSiteExpValues(const MatrixDd& firstTwoSiteOp,
                      const MatrixDd& secondTwoSiteOp, int rangeOfObservables,
                      int currentLSys, EffectiveHamiltonian& hSuperFinal,
                      std::vector<TheBlock>& leftBlocks,
                      std::vector<TheBlock>& rightBlocks,
                      std::ofstream& fileout);
void modifyHamParams(int trial = 0);

#endif
