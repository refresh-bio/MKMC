#include "parameters.h"

Params::Params()
{
	stage1ParamsTemplate.SetInputFileType(KMC::InputFileType::FASTA);

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
	if (mkmcParams.nThreads < mkmcParams.nKMCWorkers)
	{
		stage1ParamsTemplate.SetNThreads(1);
		stage2ParamsTemplate.SetNThreads(1);
	}
	else
	{
		stage1ParamsTemplate.SetNThreads(mkmcParams.nThreads / mkmcParams.nKMCWorkers);
		stage2ParamsTemplate.SetNThreads(mkmcParams.nThreads / mkmcParams.nKMCWorkers);
	}

	if (mkmcParams.maxRamGB < 2 * mkmcParams.nKMCWorkers) {
		stage1ParamsTemplate.SetMaxRamGB(2);
		stage2ParamsTemplate.SetMaxRamGB(2);
	}
	else
	{
		stage1ParamsTemplate.SetMaxRamGB(mkmcParams.maxRamGB / mkmcParams.nKMCWorkers);
		stage2ParamsTemplate.SetMaxRamGB(mkmcParams.maxRamGB / mkmcParams.nKMCWorkers);
	}
}
