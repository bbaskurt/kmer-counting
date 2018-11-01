#include "kmerBaseTypes.h"
#include "kmc.h"

/************************* The program **************************

kmerCounter is a program for counting k-mers from FASTQ files.
Input k-mers are clustered using signature which is a m-mer
(m <= k) fullfills requirements to be a signature described in
"Deorowicz, S., Kokot, M., Grabowski, S.,& Debudaj-Grabysz, A.
(2015). KMC2: fast and resource-frugal k-mer counting.
Bioinformatics, 31(10), 1569-1576." Signatures are extracted
from the sequence read from FASTQ file. Then super-kmer which
represents adjacent k-mers having same signature is written to
temporary file (bin) associated with the signature. Adjacent
k-mers can be stored in more compact owing to overlap of k-1
characters. Signature clustering provides that all occurences
of same k-mer are written to same bin even different k-mers
might have same signature. Thus, counting each bin seperately
will be enough to have total count of a k-mer that provides
huge advantage on memory usage of hash maps.

After creation of temporary files containing clustered kmers,
counter threads are created to share existing bins and count
them seperately. A counting result collector is created to merge
result to obtain top k-mers. Each counter notifies result
collector once counting operation finishes. Results are written
to "topKmers.txt" file with execution time and total number of
k-mers counted.


************************* Installation ***************************

Application does not have any 3rd party dependencies. Any
compiler with c++11 support is enough to compile source code. One
can find a simple compilation command in "compile.sh" file.


**************************** Usage *******************************

The program takes parameters three parameters such as input file,
k, and number of most frequent k-mers to be listed via command
line arguments.

./kmerCounter "filename" "kmersize" "topcount"
Example: /./kmerCounter test.fastq 30 25)

filename : fastq file name
kmersize : integer value between 1 and 256
topcount : integer value between 1 and 10000

default signature length : 6 (signatureLen = k/2, if (1 < k < 6))

******************************************************************

Author: Batuhan Baskurt
Date  : 2017-09-01

******************************************************************/

int main(int argc, char* argv[])
{
	kmerCounting::kmerParams params;

	if (argc < 4)
	{
		std::cout << "Invalid parameters" << std::endl;
		std::cout << "Usage: " << argv[0] << " file-name " << " k-mer-size " << " number-of-top-k-mers " << std::endl;
		std::cout << "Example: " << argv[0] << " big.fastq 30 25 " << std::endl;
		return -1;
	}
	else
	{
		params.filename = argv[1];
		if (params.filename.size() <= 0)
		{
			std::cout << "Invalid file path" << std::endl;
			return -1;
		}
		std::ifstream filein(params.filename);
		if (!filein.is_open())
		{
			std::cout << "Error on reading input file: " << params.filename << std::endl;
			return -1;
		}
		filein.close();

		params.kMerSize = atoi(argv[2]);
		if (params.kMerSize <= 0 || params.kMerSize > KMER_MAX_K)
		{
			std::cout << "Invalid k-mer size: " << params.kMerSize << std::endl;
			std::cout << "(k must be between 1 and " << KMER_MAX_K << ")" << std::endl;
			return -1;
		}
		if (params.kMerSize < params.signatureLen)
		{
			params.signatureLen = params.kMerSize / 2;
			if (params.signatureLen <= 0)
			{
				params.signatureLen = 1;
			}
		}
		params.topkMerCount = atoi(argv[3]);
		if (params.topkMerCount <= 0 || params.topkMerCount > KMER_MAX_TOP_COUNT)
		{
			std::cout << "Invalid top k-mer count: " << params.topkMerCount << std::endl;
			std::cout << "(Top count must be between 1 and " << KMER_MAX_TOP_COUNT << ")" << std::endl;
			return -1;
		}
	}

	kmerCounting::kmc counter(params);

	std::cin.get();

	return 0;
}