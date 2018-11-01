#include "kmerMemoryController.h"

kmerCounting::kmerMemoryController* kmerCounting::kmerMemoryController::s_instance = nullptr;
std::mutex kmerCounting::kmerMemoryController::m_instanceMutex;

kmerCounting::kmerMemoryController::kmerMemoryController()
{
}


kmerCounting::kmerMemoryController::~kmerMemoryController()
{
}

kmerCounting::kmerMemoryController* kmerCounting::kmerMemoryController::getInstance()
{
	std::lock_guard<std::mutex> lock(m_instanceMutex);

	if (s_instance == NULL)
		s_instance = new kmerMemoryController();
	return s_instance;
}


bool kmerCounting::kmerMemoryController::memoryAllocated(uint64 size)
{
	std::lock_guard<std::mutex> lock(m_instanceMutex);
	if (size < (m_maxMemory - m_totalMemoryInUse))
	{
		m_totalMemoryInUse += size;
		return true;
	}
	return false;
}

uint64 kmerCounting::kmerMemoryController::memoryFreed(uint64 size)
{
	std::lock_guard<std::mutex> lock(m_instanceMutex);
	m_totalMemoryInUse -= size;
	return m_totalMemoryInUse;
}

uint64 kmerCounting::kmerMemoryController::getFreeMemorySize()
{
	std::lock_guard<std::mutex> lock(m_instanceMutex);
	return m_maxMemory - m_totalMemoryInUse;
}