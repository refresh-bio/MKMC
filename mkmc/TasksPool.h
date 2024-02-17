#pragma once

#include <cstdint>
#include <mutex>
#include <vector>
#include <cassert>



template<typename TASK_T>
class TasksPool
{
	const std::vector<TASK_T>& tasks;

	uint32_t nextTask;
	std::mutex taskAvailableMutex;

public:
	TasksPool(const std::vector<TASK_T>& tasks) :
		tasks(tasks),
		nextTask(0)
	{}

	bool getTask(TASK_T& task);
};



template<typename TASK_T>
bool TasksPool<TASK_T>::getTask(TASK_T& task) {
	std::unique_lock<std::mutex> lck(taskAvailableMutex);
	assert(nextTask <= tasks.size());
	if (nextTask == tasks.size())
	{
		return false;
	}
	task = tasks[nextTask];
	nextTask++;
	return true;
}