#ifndef _RESULT_COLLECTOR_H
#define _RESULT_COLLECTOR_H

#include "kmerDefinitions.h"
#include <mutex>
#include <thread>
#include"kmerBaseTypes.h"

namespace kmerCounting
{
	/**
	* Singleton class for collecting kmer counting results received from counter threads.
	* Top n kmers are written to "topKmers.txt" file once counting finishes.
	*/
	class kmerResultCollector
	{
	public:
		~kmerResultCollector();
		static kmerResultCollector* getInstance();
		bool init(const kmerCounting::kmerParams& params);

		/**
		* Gets top n kmer counting result of counter thread which calls this function
		* at the end of counting operation
		*
		* @param results top n kmers of notifier thread 
		* @param totalNumberOfKmersCounted total number of kmers counted by notifier thread
		*/
		void notifyResult(const std::vector<std::pair<std::string, int64> > results, int64 totalNumberOfKmersCounted);

		/**
		* This function is called after finishing all counting operations 
		* Merges received top kmer counts to find final results. 
		*
		* @param processTime total time of whole process (in ms)
		*/
		void countingCompleted(long long processTime);

	private:
		kmerResultCollector();
		static kmerResultCollector* s_instance;
		static std::mutex m_instanceMutex;
		kmerParams m_params;
		std::vector<std::pair<std::string, int64> > m_topKmers;
		int64 m_minTop;
		int64 m_totalNumberOfKmers;
	};

}
#endif //_RESULT_COLLECTOR_H

