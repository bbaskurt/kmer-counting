#ifndef _KMER_BASE_TYPES_H
#define _KMER_BASE_TYPES_H

#include "kmerDefinitions.h"

namespace kmerCounting
{
	// Application parameters
	struct kmerParams {
		kmerParams() :
			filename(""),
			kMerSize(30),
			topkMerCount(25),
			memoryMode(false),
			signatureLen(6),		// 4^len bins for all signatures
			useCanonical(true)
			//maxX(0)
		{
			counterThreadCount = std::thread::hardware_concurrency();
			if (counterThreadCount < 1)
				counterThreadCount = 1;
		}
		std::string filename;		// input file name
		int kMerSize;				// k number
		int topkMerCount;			// number of top kmers to be listed 
		bool memoryMode;			// memory or temporary file mode (memory mode is not implemented!)
		int signatureLen;			// length of signature
		bool useCanonical;			// compares kmer and its reverse to find canonical kmer
		//int maxX;					// do not create k,x-mers if it is 0.
		int counterThreadCount;		// number of counter threads 
	};

}

#endif	// _KMER_BASE_TYPES_H
