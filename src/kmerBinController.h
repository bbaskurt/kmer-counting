#ifndef _BIN_CONTROLLER_H
#define _BIN_CONTROLLER_H

#include "kmerDefinitions.h"
#include <mutex>
#include <thread>

namespace kmerCounting
{
	/**
	* Singleton class for management of stored bins (temporary files).
	*/
	class kmerBinController
	{
	public:

		~kmerBinController();
		
		/**
		* All functions of the class are reached by using unique instance.
		*/
		static kmerBinController* getInstance();

		/**
		* Initializes instance
		*/
		bool init();

		/**
		* Increases super kmer count belonging to given bin.
		*
		* @param binId bin id
		* @param superKmerCount number of super kmers
		*/
		void addSuperKmerCountToBin(std::string binId, int64 superKmerCount);
		
		/**
		* Gets number of super kmers exist in given bin file.
		*
		* @param binId bin id
		* @return number of super kmers
		*/
		int64 getSuperKmerCount(std::string binId);

		/**
		* Increases kmer count belonging to given bin
		*
		* @param binId bin id
		* @param kmerCount number of super kmers
		*/
		void addKmerCountToBin(std::string binId, int64 kmerCount);
		
		/**
		* Gets number of kmers exist in given bin file.
		*
		* @param binId bin id
		* @return number of kmers
		*/
		int64 getKmerCount(std::string binId);

		/**
		* Increases size of given bin. This function is called after 
		* each data write to given bin file.
		*
		* @param binId bin id
		* @param binSize bin size (bytes)
		*/
		void addBinSize(std::string binId, int64 binSize);

		/**
		* Gets size of given bin file.
		*
		* @param binId bin id
		* @return bin size (bytes)
		*/
		int64 getBinSize(std::string binId);

		/**
		* Adds given file name to list of given bin id.
		*
		* @param binId bin id
		* @param filename file name in which super kmers of given bin id is stored.
		* @return number of files exist for the bin id.
		*/
		int addFileToBin(std::string binId, const std::string& filename);

		/**
		* Gets related file list of given bin id.
		*
		* @param binId bin id
		* @return file list, if exists.
		*		  empty list, otherwise.
		*/
		std::vector<std::string> getFileListOfBin(std::string binId);

		/**
		* Fills all existing bin ids except supplementary ones(having '_' in the id) to given list 
		*
		* @param binIdList bin id list
		*/
		void getAllBinIds(std::vector<int32>& binIdList);

		/**
		* Gets total number of kmers read from input file
		*
		* @return total number of kmers read from input file
		*/
		int64 getTotalNumberOfKmers();

		/**
		* Gets maximum number of kmers storing in a single bin.
		*
		* @return kmer count
		*/
		int64 getMaxKmerCountForBin();

		/**
		* Removes all temporary files.
		*/
		void removeTempFiles();

	private:
		kmerBinController();
		static kmerBinController* s_instance;
		static std::mutex m_instanceMutex;
		// map to store file sizes (format: <bin id, data size>)
		// bin id is stored in string considering multiple files for large bins (i.e 0097, 0097_1, 0097_2, etc.)
		std::unordered_map<std::string, int64> m_binSizeMap;
		// map to store filenames of bins bigger than maximum single file size (format: <bin id, data size>)
		// Stores only supplementary files.
		// (Normally, all bins are stored in file names with numeric representation of their signature.)
		std::unordered_map<std::string, std::vector<std::string> > m_binsWithMultipleFiles;
		// map to store kmer count of each bin (format: <bin id, number of kmers>)
		std::unordered_map<std::string, int64> m_kmerCountForBin;
		// map to store super kmer count of each bin (format: <bin id, number of super kmers>)
		std::unordered_map<std::string, int64> m_superkmerCountForBin;
		int64 m_maxKmerCountForBin;
	};

}
#endif //_BIN_CONTROLLER_H

