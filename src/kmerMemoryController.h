#ifndef _MEMORY_CONTROLLER_H
#define _MEMORY_CONTROLLER_H

#include "kmerDefinitions.h"
#include <mutex>
#include <thread>

namespace kmerCounting
{
	/**
	* Singleton class for controlling memory currently being used by the application
	*/
	class kmerMemoryController
	{
	public:
		~kmerMemoryController();
		static kmerMemoryController* getInstance();
		bool init();

		/**
		* Checks if there is free memory of given size to allocate. 
		* Increases memory in use if it is available, otherwise returns false.
		* This function is called before memory allocation.
		*
		* @param size number of bytes to allocate.
		* @return true, if given size is available
		*		  false, otherwise
		*/
		bool memoryAllocated(uint64 size);

		/**
		* Updates memory in use according to given size of memory freed recently
		*
		* @param size number of bytes that has been freed.
		* @return total number of bytes currently in use
		*/
		uint64 memoryFreed(uint64 size);

		/**
		* Gets number of bytes available that application is allowed to use.
		*
		* @return total number of available bytes.
		*/
		uint64 getFreeMemorySize();

	private:
		kmerMemoryController();
		static kmerMemoryController* s_instance;
		static std::mutex m_instanceMutex;
		uint64 m_totalMemoryInUse;		// bytes
		uint64 m_maxMemory;				// max number of memory that application is allowed to use (bytes)

	};

}
#endif //_MEMORY_CONTROLLER_H

