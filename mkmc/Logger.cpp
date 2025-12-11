#include "Logger.h"
#include <iostream>



Logger& Logger::Inst() {
	static Logger inst;
	return inst;
}


bool Logger::Enable(int _cerrVerbosityLevel, const std::string& logFileName) {
	logStream.open(logFileName);
	if (!logStream.is_open())
		std::cerr << "Warning: cannot create log file " << logFileName << ". Log will not be safed to the file." << std::endl;

	is_enabled = true;
	cerrVerbosityLevel = _cerrVerbosityLevel;

	return true;
}


void Logger::Disable() {
	is_enabled = false;
	logStream.close();
}


void Logger::Log(const std::string& msg, int msgVerbosityLevel) {
	if (!is_enabled)
		return;

	std::lock_guard<std::mutex> lck(mtx);
	if (msgVerbosityLevel <= cerrVerbosityLevel)
		std::cerr << msg << "\n";

	if (logStream.is_open())
		logStream << msg << std::endl; // Flush immediately to get information before MKMC finish/filling the buffer
}
