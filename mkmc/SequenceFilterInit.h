#pragma once

#include <vector>
#include <fstream>
#include <sstream>
#include <filesystem>
#include "parameters.h"
#include "time.hpp"



class SequenceFilterInit
{
	const Params& params;
	Timer& sequence_filter_init;

	bool getNotEmptyLine(std::istream& stream, std::string& outLine);

	bool isFastaOrMultiFasta();
	void convertTxtToFasta();

public:
	SequenceFilterInit(const Params& params, Timer& sequence_filter_init) :
		params(params),
		sequence_filter_init(sequence_filter_init)
	{}

	bool prepareKmersSequencesToFilter();
};
