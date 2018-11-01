#ifndef _KMC_H
#define _KMC_H

#include "kmerBaseTypes.h"

namespace kmerCounting
{
	class kmerExtractor;
	class kmerCounter;
	class kmc
	{
	public:
		/**
		* Starts kmer counting operation. First, kmerExtractor which extracts superkmers and writes them to temporary files 
		* is created. Then counter threads are created according to number of computer's core. Total number of temporary 
		* files are shared to counter threads. Each thread has own kmerCounter object which counts corresponding kmers and 
		* notifies kmerResultCollector once finishing its operation. kmerResultCollector merges all results received from
		* counter threads to obtain final counts. Top kmers are written in both console output and "topKmers.txt" file in current 
		* directory.
		*
		* @param params algorithm parameters determined by user.
		*/
		kmc(const kmerCounting::kmerParams& params);

		~kmc();

	private:		
		kmerParams m_params;
		kmerExtractor *m_kmerExtractor;
	};
}

#endif	// _KMC_H


