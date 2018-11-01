#ifndef _KMER_EXTRACTOR_H
#define _KMER_EXTRACTOR_H

#include "kmerDefinitions.h"
#include "kmerBaseTypes.h"

namespace kmerCounting
{
	class fastqReader;
	class kmerExtractor
	{
	public:

		/**
		* k-mers of input file are extracted and stored in temporary files using kmerExtractor.
		* k-mers are clustered according to their signatures. A temporary file containing kmers 
		* of each signature is created. Extracted k-mers are first stored in a map in memory until
		* reaching a pre-defined size (KMER_BIN_BUFFER_SIZE) to write temporary file.
		*
		* @param params algorithm parameters determined by user.
		*/
		kmerExtractor(const kmerCounting::kmerParams& params);

		~kmerExtractor();
		
		/**
		* Extract k-mers of given file name in params.
		*/
		void extractSuperkmers();

		/**
		* Puts given super k-mer to related bin. 
		*
		* @param binNo bin id
		* @param seq data pointer to super k-mer
		* @return n length of super k-mer
		*/
		void putExtendedKmer(uint32 binNo, char* seq, uint32 n);

	private:

		/**
		* Gets temporary file name of given bin id. Current size of temporary file and 
		* number of bytes to be written are considered to make decision of creating 
		* supplementary file for the bin.
		*
		* @param binNo bin id
		* @param bytes bytes to be written
		* @return file name
		*/
		std::string getFileName(uint32 binNo, int32 bytes);

		/**
		* Increases number of super k-mers storing in local map belonging to given bin. 
		*
		* @param binNo bin id
		* @param increment amount of increment
		* @return total number
		*/
		uint32 increaseSuperKmerCount(uint32 binNo, uint32 increment);

		/**
		* Increases number of k-mers storing in local map belonging to given bin.
		*
		* @param binNo bin id
		* @param increment amount of increment
		* @return total number
		*/
		uint32 increaseKmerCount(uint32 binNo, uint32 increment);

	private:
		kmerParams m_params;
		// stores each possible signature
		int32* m_signatureMap;				// use a seperate file (bin) for each kmerSignature
		// number of possible signatures according to given signature length
		int32 m_mapSize;
		// map to store super kmers (format: <bin id, <data pointer, data size> >)
		std::unordered_map<int32, std::pair<uchar*, int32> > m_bins;
		// map to store number of super kmers for current bin buffer(format: <bin id, number of super kmers >)
		// value of each buffer is set to zero on each data transfer from buffer to file
		std::unordered_map<int32, uint32> m_superKmersCount;
		// map to store number of kmers for current bin buffer (format: <bin id, number of kmers >)
		// value of each buffer is set to zero on each data transfer from buffer to file
		std::unordered_map<int32, uint32> m_kmerCount;
		//uint32 m_superKmersCount;
		//uint32 m_recsCount;
		uint32 m_plusXRecsCount;
		//uint32 m_maxX;

		fastqReader* m_fastqReader;
	};
}
#endif	// _KMER_EXTRACTOR_H
