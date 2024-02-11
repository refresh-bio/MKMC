#pragma once

#include <cstdint>
#include <mutex>

class TasksPool
{
	uint32_t nTasks;
	uint32_t nextTask;
	std::mutex taskAvailableMutex;

public:
	TasksPool(const uint32_t nTasks) :
		nTasks(nTasks),
		nextTask(0)
	{}

	bool getTask(uint32_t& task);
};
