#ifndef _KMER_UTILS_H
#define _KMER_UTILS_H

#include "kmerDefinitions.h"

namespace kmerCounting
{
	class kmerUtils
	{
	public:

		/**
		* Commonly used kmer functions.
		*/
		kmerUtils();

		~kmerUtils();

		/**
		* Checks if given mmer satisfy requirement of being a kmerSignature.
		*
		* @param mmer mmer
		* @param len mmer length
		* @return result
		*/
		static bool isAllowed(uint32 mmer, uint32 len);

		/**
		* Gets reverse complement of given mmer
		*
		* @param mmer mmer
		* @param len mmer length
		* @return reverse complement of given mmer
		*/
		static uint32 getRev(uint32 mmer, uint32 len);

		/**
		* Gets symbol (A,C,G,T) of given numberic value
		*
		* @param value numeric value
		* @return char correspondence of value
		*/
		static char getSymbol(char value);
	};
}
#endif	// _KMER_UTILS_H

