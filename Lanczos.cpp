#include "d.h"
#include "main.h"

extern "C"
{
    void dstemr_(char* JOBZ, char* RANGE, int* N, double* D, double* E,
                 double* VL, double* VU, int* IL, int* IU, int* M, double* W,
                 double* Z, int* LDZ, int* NZC, int* ISUPPZ, bool* TRYRAC,
                 double* WORK, int* LWORK, int* IWORK, int* LIWORK, int* INFO);
};

using namespace Eigen;

std::pair<VectorXd, double> lanczos(const MatrixXd& mat, double tolerance)
{
    std::vector<double> a,
                        b;
    int n = mat.rows();
    VectorXd x = VectorXd::Random(n);
    MatrixXd basisVecs = x / x.norm();              // initial seed
    x.noalias() = mat * basisVecs;
    a.push_back(basisVecs.col(0).dot(x));
    b.push_back(0.);
    VectorXd oldGS,
             newGS = basisVecs;
    int i = 0;                                      // iteration counter
    char JOBZ = 'V',                                // define dstemr arguments
         RANGE = 'I';
    int N = 1;
    std::vector<double> D,
                        E;
    double VL,
           VU;
    int IL = 1,
        IU = 1,
        M;
    std::vector<double> W;
    VectorXd Z;
    int LDZ,
        NZC = 1;
    std::vector<int> ISUPPZ;
    ISUPPZ.reserve(2);
    bool TRYRAC = false;
    double optLWORK;
    std::vector<double> WORK;
    int LWORK,
        optLIWORK;
    std::vector<int> IWORK;
    int LIWORK,
        INFO;
    do
    {
        i++;
        oldGS = newGS;
        
        // Lanczos stage 1: Lanczos iteration
        x -= a[i - 1] * basisVecs.col(i - 1);
        b.push_back(x.norm());
        basisVecs.conservativeResize(NoChange, i + 1);
        basisVecs.col(i) = x / b[i];
        x.noalias() = mat * basisVecs.col(i) - b[i] * basisVecs.col(i - 1);
        a.push_back(basisVecs.col(i).dot(x));
        b.reserve(i + 2);                  // allocate an extra slot for dstemr
        
        // Lanczos stage 2: diagonalize tridiagonal matrix
        N++;
        D = a;
//        E.assign(b.begin() + 1, b.end());
//        E = std::vector<double>(b.begin() + 1, b.end());
        E = b;
        E.erase(E.begin());
        W.reserve(N);
        Z.resize(N);
        LDZ = N;
        LWORK = -1;
        LIWORK = -1;
        dstemr_(&JOBZ, &RANGE, &N, D.data(), E.data(), &VL, &VU, &IL, &IU, &M,
                W.data(), Z.data(), &LDZ, &NZC, ISUPPZ.data(), &TRYRAC,
                &optLWORK, &LWORK, &optLIWORK, &LIWORK, &INFO);
                                    // query for optimal workspace allocations
        LWORK = int(optLWORK);
        WORK.reserve(LWORK);
        LIWORK = optLIWORK;
        IWORK.reserve(LIWORK);
        dstemr_(&JOBZ, &RANGE, &N, D.data(), E.data(), &VL, &VU, &IL, &IU, &M,
                W.data(), Z.data(), &LDZ, &NZC, ISUPPZ.data(), &TRYRAC,
                WORK.data(), &LWORK, IWORK.data(), &LIWORK, &INFO);
        newGS = basisVecs * Z;
    } while(std::min((newGS - oldGS).norm(), (newGS + oldGS).norm()) > tolerance
            && N <= n);
    if(N > n)
    {
        std::cerr << "Lanczos algorithm failed to converge." << std::endl;
        exit(EXIT_FAILURE);
    };
    return std::pair<VectorXd, double>(newGS, W.front());
};
