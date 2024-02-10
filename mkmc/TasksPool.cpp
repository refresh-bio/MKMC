#include "TasksPool.h"
#include <cassert>

bool TasksPool::getTask(uint32_t& task) {
	std::unique_lock<std::mutex> lck(taskAvailableMutex);
	assert(nextTask <= nTasks);
	if (nextTask == nTasks)
	{
		return false;
	}
	task = nextTask++;
	return true;
}
