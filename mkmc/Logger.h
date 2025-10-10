#pragma once
#include <mutex>
#include <string>
#include <fstream>


//thread safe singleton
class Logger {

	mutable std::mutex mtx;
	bool is_enabled = false;

	int cerrVerbosityLevel = 0;

	std::ofstream logStream;

public:
	static Logger& Inst();
	bool Enable(int _cerrVerbosityLevel, const std::string& logFileName);
	void Disable();
	void Log(const std::string& msg, int msgVerbosityLevel = 0);
};
