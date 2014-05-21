#ifndef FINALSUPERBLOCK_H
#define FINALSUPERBLOCK_H

#include <map>
#include "TheBlock.h"

typedef std::vector<std::pair<obsMatrixD_t, int>,
                    Eigen::aligned_allocator<std::pair<obsMatrixD_t, int>>>
                    opsVec;
typedef std::map<int, obsMatrixD_t, std::less<int>,
                 Eigen::aligned_allocator<std::pair<const int, obsMatrixD_t>>>
                 opsMap;

class FinalSuperblock
{
    public:
        double gsEnergy;                                 // ground-state energy
        
        FinalSuperblock(const MatrixX_t& matFinal,
                        const std::vector<int>& qNumList,
                        const std::vector<int>& compQNumList,
                        const stepData& data,
                        const rmMatrixX_t& psiGroundIn, int mSFinal,
                        int mEFinal, int skips);
        double expValue(const opsVec& ops, std::vector<TheBlock>& leftBlocks,
                        std::vector<TheBlock>& rightBlocks);
        // calculates exectation value of a combination of single-site operators
    
    private:
        int lSupFinal;                                     // final system size
        rmMatrixX_t psiGround;                 // final superblock ground state
        int lSFinal,                            // final length of system block
            lEFinal,                       // final length of environment block
            mSFinal,           // final number of states stored in system block
            mEFinal,      // final number of states stored in environment block
            skips;                // number of edge sites in the position basis
        
        void placeOp(const std::pair<obsMatrixD_t, int>& op, opsMap& blockSide,
                     bool systemSide);
                    // assign each one-site observable to the appropriate block
        obsMatrixX_t rhoBasisRep(const opsMap& blockOps,
                                 std::vector<TheBlock>& blocks, int blockSize)
                                 const;
                  // converts single-site operators into the system block basis
};

#endif
