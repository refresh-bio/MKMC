#include "parameters.h"
#include <iostream>

Params::Params()
{
	stage1Params.SetInputFileType(KMC::InputFileType::FASTQ);
	stage1Params.SetNBins(64);

	stage2Params.SetCutoffMin(1);
	stage2Params.SetCutoffMax(static_cast<uint64_t>(4E9));
	stage2Params.SetCounterMax(65535);

	static KMC::NullPercentProgressObserver nullPercentProgressObserver;
	static KMC::NullProgressObserver nullProgressObserver;
	stage1Params.SetPercentProgressObserver(&nullPercentProgressObserver);
	stage1Params.SetProgressObserver(&nullProgressObserver);
}

void Params::setKMCParams()
{
	bool mKMCWorkersReduced = false;
	if (mkmcParams.nKMCWorkers > mkmcParams.nThreads - 1)
	{
		mkmcParams.nKMCWorkers = mkmcParams.nThreads - 1;
		mKMCWorkersReduced = true;
	}

	stage1Params.SetNThreads(mkmcParams.nThreads / mkmcParams.nKMCWorkers);
	stage2Params.SetNThreads(mkmcParams.nThreads / mkmcParams.nKMCWorkers);

	if (mkmcParams.nKMCWorkers * 2 > mkmcParams.maxRamGB)
	{
		mkmcParams.nKMCWorkers = mkmcParams.maxRamGB / 2;
		mKMCWorkersReduced = true;
	}

	if (mkmcParams.nKMCWorkersUserSet && mKMCWorkersReduced)
	{
		std::cerr << "Warning: number of workers is too huge, reduced to " << mkmcParams.nKMCWorkers << std::endl;
	}

	stage1Params.SetMaxRamGB(mkmcParams.maxRamGB / mkmcParams.nKMCWorkers);
	stage2Params.SetMaxRamGB(mkmcParams.maxRamGB / mkmcParams.nKMCWorkers);
}
