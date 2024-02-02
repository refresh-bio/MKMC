#include "parameters.h"
#include <iostream>

Params::Params()
{
	stage1ParamsTemplate.SetInputFileType(KMC::InputFileType::FASTQ);

	stage2ParamsTemplate.SetCutoffMin(1);
	stage2ParamsTemplate.SetCutoffMax(static_cast<uint64_t>(4E9));
	stage2ParamsTemplate.SetCounterMax(65535);

	static KMC::NullPercentProgressObserver nullPercentProgressObserver;
	static KMC::NullProgressObserver nullProgressObserver;
	stage1ParamsTemplate.SetPercentProgressObserver(&nullPercentProgressObserver);
	stage1ParamsTemplate.SetProgressObserver(&nullProgressObserver);
}

void Params::setKMCParams()
{
	bool mKMCWorkersReduced = false;
	if (mkmcParams.nKMCWorkers > mkmcParams.nThreads - 1)
	{
		mkmcParams.nKMCWorkers = mkmcParams.nThreads - 1;
		mKMCWorkersReduced = true;
	}

	stage1ParamsTemplate.SetNThreads(mkmcParams.nThreads / mkmcParams.nKMCWorkers);
	stage2ParamsTemplate.SetNThreads(mkmcParams.nThreads / mkmcParams.nKMCWorkers);

	if (mkmcParams.nKMCWorkers * 2 > mkmcParams.maxRamGB)
	{
		mkmcParams.nKMCWorkers = mkmcParams.maxRamGB / 2;
		mKMCWorkersReduced = true;
	}

	if (mkmcParams.nKMCWorkersUserSet && mKMCWorkersReduced)
	{
		std::cerr << "Warning: number of workers is too huge, reduced to " << mkmcParams.nKMCWorkers << std::endl;
	}

	stage1ParamsTemplate.SetMaxRamGB(mkmcParams.maxRamGB / mkmcParams.nKMCWorkers);
	stage2ParamsTemplate.SetMaxRamGB(mkmcParams.maxRamGB / mkmcParams.nKMCWorkers);

	stage1ParamsTemplate.SetRamOnlyMode(true);
}
