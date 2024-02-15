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


class PercentProgress
{
	size_t curr_val;
	size_t max_val;	
	int32_t curr_percent;
	bool show_progress;
public:
	PercentProgress(size_t max_val, bool show_progress)
		:
		curr_val(0),
		max_val(max_val),
		curr_percent(-1),
		show_progress(show_progress)
	{
		
	}
	
	void NotifyProgress(uint64_t val)
	{
		if (!show_progress)
			return;

		curr_val += val;
		int32_t new_percent = 0;
		if (max_val)
			new_percent = static_cast<int32>((curr_val * 100) / max_val);
		else
			new_percent = 100;
		if (new_percent > curr_percent)
		{
			curr_percent = new_percent;
			std::cerr << "\r" << curr_percent << "%";
			if (new_percent == 100)
				std::cerr << "\n";
		}
	}
};
