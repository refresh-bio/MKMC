#pragma once

#include "parameters.h"

class Finish
{
	const Params& params;
	
public:
	Finish(const Params& params) :
		params(params)
	{}

	void finishProcessing();
};