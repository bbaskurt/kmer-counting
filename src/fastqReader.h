#ifndef _FASTQ_READER_H
#define _FASTQ_READER_H

#include "kmerDefinitions.h"
#include <stdio.h>
#include <string>
#include <iostream>

namespace kmerCounting
{
	class fastqReader
	{
	public:

		/**
		* I/O operations of fastq file are operated in this class.
		*
		* @param memoryMode data operations are done on memory instead of file, if it is true.
		*/
		fastqReader(bool memoryMode = false);

		~fastqReader();

		/**
		* Gets a single sequence from input file. Input file is read in partitions of KMER_FILE_READER_PART_SIZE.
		* So sequence is first being tried to taken from memory. If there is no enough data remained in memory, 
		* new partition is read from the file.
		*
		* @param seq input sequence
		* @param seqSize sequence size (bytes)
		* @return status
		*/
		bool getSequence(char *seq, uint32 &seqSize);

		/**
		* Opens given file.
		*
		* @param fName file name
		*/
		void open(const std::string& fName);

		/**
		* Closes current file.
		*
		* @return operation result
		*/
		int close();

		/**
		* Reads binary data from input file according to given data size and count.
		*
		* @param ptr data pointer to be filled
		* @param size element size (bytes)
		* @param count element count
		* @return size of read data
		*/
		size_t read(uchar * ptr, size_t size, size_t count);

		/**
		* Writes binary data to current file according to given data size and count.
		*
		* @param ptr data pointer to be written
		* @param size element size (bytes)
		* @param count element count
		* @return size of written data
		*/
		size_t write(const uchar * ptr, size_t size, size_t count);

	private:

		/**
		* Checks if part on the memory has enough characters to read. 
		* Otherwise reads more data from input file.
		*
		* @return operation result
		*/
		bool isPartAvailable();

	private:
		bool m_memoryMode;
		FILE* m_file;
		int64 m_filePos;
		int64 m_fileSize;
		int64 m_partPos;
		int64 m_partSize;
		uchar* m_part;
		char m_codes[MAX_NUMBER_OF_SYMBOLS];
	};
}
#endif //_FASTQ_READER_H


