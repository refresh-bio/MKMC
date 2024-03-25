#pragma once
#include <mutex>
#include <string>
#include <iostream>


//thread safe singleton
class Logger {

	mutable std::mutex mtx;
	bool is_enabled = false;

public:
	static Logger& Inst() {
		static Logger inst;
		return inst;
	}
	void Enable() {
		is_enabled = true;
	}
	void Disable() {
		is_enabled = false;
	}
	void Log(const std::string& msg) const {
		if (!is_enabled)
			return;
		std::lock_guard<std::mutex> lck(mtx);
		std::cerr << msg << "\n";
	}
};
