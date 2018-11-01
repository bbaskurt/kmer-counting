#ifndef _KMER_SIGNATURE_H
#define _KMER_SIGNATURE_H

#include "kmerDefinitions.h"

namespace kmerCounting
{

	class kmerSignature
	{
	public:

		/**
		* kmer signature which is mmer (m<=k) satisfying requirements described in "Deorowicz, S., Kokot, M., Grabowski, S.,
		* & Debudaj-Grabysz, A. (2015). KMC 2: fast and resource-frugal k-mer counting. Bioinformatics, 31(10), 1569-1576."
		* Signature is stored by its numeric representation in which each symbol (A,C,G,T) is stored in 2 bits.
		*
		* @param signatureLen signature length determined by user.
		*/
		kmerSignature(const int& signatureLen);

		~kmerSignature();
		
		/**
		* Inserts given symbol to signature.
		*
		* @param seq symbol
		*/
		void insert(char* seq);

		/**
		* Inserts given symbol to signature.
		*
		* @param seq symbol
		*/
		void insert(uchar symb);

		/**
		* Gets numeric representation of the signature
		*
		* @return signature
		*/
		uint32 get() const;

		/**
		* Clears signature
		*/
		void clear();
		
		/**
		* Sets signature
		*
		* @param sig signature
		*/
		void set(const kmerSignature& sig);

		bool operator==(const kmerSignature& x);
		bool operator<(const kmerSignature& x);
		bool operator<=(const kmerSignature& x);

	private:

		/**
		* Creates a look-up table containing all possible signatures corresponding to given signature length.
		* Value of each signature is numerically minimum value of its reverse complement and itself. 
		* If the signature is not valid according to rules described in the referenced article, a special value 
		* is assigned. Same special value is used for all invalid signatures that means all kmers do not have 
		* valid signature will be written to same bin by kmerExtractor.
		*
		* @param sig signature
		*/
		void init_norm(uint32* norm, uint32 len);

	private:
		uint32 m_len;
		uint32* m_norm;
		uint32 m_currentVal;
		uint32 m_str;
		uint32 m_mask;

		static uint32 norm1[1 << 2];
		static uint32 norm2[1 << 4];
		static uint32 norm3[1 << 6];
		static uint32 norm4[1 << 8];
		static uint32 norm5[1 << 10];
		static uint32 norm6[1 << 12];
		static uint32 norm7[1 << 14];
		static uint32 norm8[1 << 16];
	};

}
#endif	// _KMER_SIGNATURE_H
