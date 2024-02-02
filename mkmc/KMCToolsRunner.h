#pragma once

#if defined(WIN32) || defined(_WIN32)
#include <Windows.h>
#else
#include <sys/stat.h>
#endif

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <string>
#include "parameters.h"
#include "KMCRunner.h"



#if defined(WIN32) || defined(_WIN32) // Windows
#define KMC_TOOLS_EXECUTABLE_NAME "kmc_tools.exe"
#else // Linux
#define KMC_TOOLS_EXECUTABLE_NAME "kmc_tools"
#endif

class KMCToolsRunner {
	const Params& params;

	void operator()(TasksPool& tasksPool);

	bool checkToolsRequired(const std::string& kmcOutputFile);

public:
	KMCToolsRunner(const Params& params) :
		params(params)
	{}

	void runKMCToolsParallel();

#if defined(WIN32) || defined(_WIN32) // Windows

private:
	// buffer - output for child process' stdout and stderr, if NULL, the result is sent to stdout
	bool runCommand(std::string command, const std::string& args, unsigned long& processResult, char* buffer = NULL, const int bufferSize = 0);

	bool createPipe(HANDLE& pipeReadHandle, HANDLE& pipeWriteHandle);
	void destroyPipe(HANDLE& pipeReadHandle, HANDLE& pipeWriteHandle);

#else // Linux

private:
	// buffer - output for child process' stdout and stderr, if NULL, the result is sent to stdout
	bool runCommand(std::string command, const std::string& args, unsigned long& processResult, char* buffer = NULL, const int bufferSize = 0);
#endif
};