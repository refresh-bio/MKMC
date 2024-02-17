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
#include "TasksPool.h"



#define RUN_SINGLE_TOOLS



#if defined(WIN32) || defined(_WIN32) // Windows
#define KMC_TOOLS_EXECUTABLE_NAME "kmc_tools.exe"
#else // Linux
#define KMC_TOOLS_EXECUTABLE_NAME "kmc_tools"
#endif

class KMCToolsRunner {
	struct TaskData
	{
		std::string inputFile;
		std::string outputFile;
	};
	std::vector<TaskData> tasksData;

	const Params& params;

	TasksPool<TaskData> tasksPool;
	void operator()();

	bool checkToolsRequired(const std::string& kmcOutputFile);

public:
	KMCToolsRunner(const Params& params) :
		params(params),
		tasksPool(tasksData)
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