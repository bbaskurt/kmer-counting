#ifndef _KMER_COUNTER_H
#define _KMER_COUNTER_H

#include "kmerBaseTypes.h"

namespace kmerCounting
{
	class kmer;
	class kmerCounter
	{
	public:

		/**
		* Processes temporary files containing super kmers to extract and count kmers. Counting results are 
		* stored in local map. Top counts are notified to kmerResultCollector, once counting finishes.
		*
		* @param params application parameters determined by user.
		*/
		kmerCounter(const kmerCounting::kmerParams& params);

		~kmerCounter();
		
		/**
		* Processes given bins.
		*
		* @param binIds id(numeric representation of related signature) list of bins to be processed. 
		*/
		void processBins(const std::vector<int32> binIds);

	private:

		/**
		* Processes given bin.
		*
		* @param binNo id(numeric representation of related signature) of bin to be processed.
		*/
		void processBin(const int32& binId);

		/**
		* Increases given kmer's occurance count
		*
		* @param kmer alphabetical representation of kmer
		*/
		void countKmer(const std::string& kmer);

		/**
		* Updates top kmers' count with current kmerCounts
		*/
		void updateTopKmers();

	private:
		kmerParams m_params;
		int64 m_maxKmerSizeForBin;
		// <kmer string, occurance count>
		std::unordered_map<std::string, int64> m_kmerCounts;
		std::vector<std::pair<std::string, int64> > m_topKmers;
		int64 m_minTop;
	};
}

#endif	// _KMER_COUNTER_H

