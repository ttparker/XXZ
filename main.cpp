#include <time.h>
#include <fstream>
#include "FreeFunctions.h"
#include "ESolver.h"

using namespace Eigen;

int main()
{
    clock_t start = clock();

	// **************** begin modifiable parameters
	const int numberOfTrials = 1;
	int lSys = 20;								// system length - must be even
	std::vector<double> couplingConstants;
	couplingConstants.reserve(2);
	couplingConstants.push_back(5.);				// J_{xy}
	couplingConstants.push_back(2.);				// J_z
	int targetQNum = lSys * 1 / 5;				// targeted fractional average
				// magnetization per site - product with lSys must be integer
	int mMax = 16,								// max number of stored states
		nSweeps = 2;						// number of sweeps to be performed
    double groundStateErrorTolerance = 1e-6;

	bool calcOneSiteExpValues = true, // calculate single-site expectation values?
		 calcTwoSiteExpValues = true; // calculate two-site expectation values?
	int rangeOfObservables = 20;	// number of sites in center at which to
			// calculate observables (must be even and less than currentLSys)
	MatrixDd oneSiteOp,
			 firstTwoSiteOp,
			 secondTwoSiteOp;
	MatrixDd sigmax;
	sigmax << 0., 1.,
			  1., 0.;
	MatrixDd sigmaz;
	sigmaz << 1., 0.,
			  0., -1.;
	oneSiteOp = sigmaz,
	firstTwoSiteOp = sigmax,
	secondTwoSiteOp = sigmax;
	// **************** end modifiable parameters

	Hamiltonian ham;		    		// initialize the system's Hamiltonian
	std::ofstream fileout;
	fileout.open("Output/Output", std::ios::out);
	if (!fileout)
	{
		std::cout << "Couldn't open output file." << std::endl;
		exit(1);
	};
	for(int trial = 0; trial < numberOfTrials; trial++)
	{
        ham.setParams(lSys, couplingConstants, targetQNum);
		std::cout << "Trial " << trial << ":" <<std::endl;
		fileout << "Trial " << trial << ":" <<std::endl;
        int skips = 0;
        for(int runningKeptStates = d * d; runningKeptStates <= mMax; skips++)
            runningKeptStates *= d; // find how many edge sites can be skipped
		std::vector<TheBlock> blocks(ham.lSys - 3 - skips);	// initialize system
		blocks[0] = TheBlock(ham, mMax);	// initialize the one-site block
        std::cout << "Performing iDMRG..." << std::endl;
        for(int site = 0; site < skips; site++)
            blocks[site + 1] = blocks[site].nextBlock(ham);       // initial ED
        Sector::lancTolerance = groundStateErrorTolerance
                                * groundStateErrorTolerance / 2;
        int lSFinal = ham.lSys / 2 - 1;     // final length of the system block
        for(int site = skips, end = lSFinal - 1; site < end; site++)
            blocks[site + 1] = blocks[site].nextBlock(ham, false, true, site);
        if(nSweeps != 0)
            std::cout << "Performing fDMRG..." << std::endl;
        for(int i = 1; i <= nSweeps; i++)           // perform the fDMRG sweeps
        {
            for(int site = lSFinal - 1, end = ham.lSys - 4 - skips; site < end; site++)
                blocks[site + 1] = blocks[site].nextBlock(ham, false, false, site,
                                                          blocks[ham.lSys - 4 - site],
                                                          blocks[ham.lSys - 5 - site]);
            blocks[skips].reflectPredictedPsi();
                               // reflect the system to reverse sweep direction
            for(int site = skips, end = lSFinal - 1; site < end; site++)
                blocks[site + 1] = blocks[site].nextBlock(ham, false, false, site,
                                                          blocks[ham.lSys - 4 - site],
                                                          blocks[ham.lSys - 5 - site]);
            std::cout << "Sweep " << i << " complete." << std::endl;
        };
        
		EffectiveHamiltonian hSuperFinal = blocks[lSFinal - 1]
                                 .createHSuperFinal(ham, skips);
											// calculate ground-state energy
		fileout << "Ground state energy density = "
				<< hSuperFinal.gsEnergy / ham.lSys << std::endl	<< std::endl;
		std::cout << "Calculating observables..." << std::endl;
		if(calcOneSiteExpValues)	// calculate one-site expectation values
			oneSiteExpValues(oneSiteOp, rangeOfObservables, ham.lSys,
							 hSuperFinal, blocks, fileout);
		if(calcTwoSiteExpValues)	// calculate two-site expectation values
			twoSiteExpValues(firstTwoSiteOp, secondTwoSiteOp,
							 rangeOfObservables, ham.lSys, hSuperFinal, blocks,
							 fileout);
		std::cout << std::endl;
		fileout << std::endl;
	};

	clock_t stop = clock();
    std::cout << "Done." << std::endl;
	fileout << "Elapsed time (s): " << (double)(stop - start)/CLOCKS_PER_SEC
			<< std::endl;
	fileout.close();

	return 0;
}
