#ifndef ESOLVER_H
#define ESOLVER_H

class HamSolver
{
    public:
        rmMatrixX_t lowestEvec;
        double lowestEval;
        
        HamSolver(const MatrixX_t& mat, const std::vector<int>& hSprimeQNumList,
                  const std::vector<int>& hEprimeQNumList, int targetQNum,
                  rmMatrixX_t& bigSeed, double lancTolerance);
    
    private:
        std::vector<int> hSprimeQNumList,
                         hEprimeQNumList;
        int targetQNum;
        
        double lanczos(const MatrixX_t& mat, rmMatrixX_t& seed,
                       double lancTolerance);
     // changes input seed to ground eigenvector - make sure seed is normalized
    
    friend class DMSolver;
};

class DMSector
{
    public:
        DMSector() {};
    
    private:
        Eigen::SelfAdjointEigenSolver<MatrixX_t> solver;      // DM eigensystem
        MatrixX_t sectorMat;                                 // sector operator
        std::vector<int> positions;
                              // which rows and columns of matrix are in sector
        int multiplicity;                       // size of this symmetry sector
        int sectorColumnCounter;           // tracks which sector eigenvector
                                           // to fill into a matrix eigenvector
        
        DMSector(const MatrixX_t& sectorMat, const std::vector<int>& positions);
        void solveForAll();
        VectorX_t nextHighestEvec(int fullMatrixSize);
        VectorX_t filledOutEvec(VectorX_t sectorEvec, int fullMatrixSize) const;
    
    friend class DMSolver;
};

class DMSolver
{
    private:
        MatrixX_t highestEvecs;
        std::vector<int> highestEvecQNums;
        double truncationError;
        
        DMSolver(const HamSolver hSuperSolver, int maxEvecsToKeep);
    
    friend class TheBlock;
};

#endif
