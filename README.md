k-mer Counting
=
kmerCounter is a program for counting k-mers from FASTQ files. Input k-mers are clustered using signature which is a m-mer (m <= k) fullfills requirements to be a signature described in "Deorowicz, S., Kokot, M., Grabowski, S.,& Debudaj-Grabysz, A. (2015). KMC2: fast and resource-frugal k-mer counting. Bioinformatics, 31(10), 1569-1576." Signatures are extracted from the sequence read from FASTQ file. Then super-kmer which represents adjacent k-mers having same signature is written to temporary file (bin) associated with the signature. Adjacent k-mers can be stored in more compact owing to overlap of k-1 characters. Signature clustering provides that all occurences of same k-mer are written to same bin even different k-mers might have same signature. Thus, counting each bin seperately will be enough to have total count of a k-mer that provides huge advantage on memory usage of hash maps.

After creation of temporary files containing clustered kmers, counter threads are created to share existing bins and count them seperately. A counting result collector is created to merge result to obtain top k-mers. Each counter notifies result collector once counting operation finishes. Results are written to "topKmers.txt" file with execution time and total number of k-mers counted. 

Installation 
=
Application does not have any 3rd party dependencies. Any compiler with c++11 support is enough to compile source code. One can find a simple compilation command in "compile.sh" file. 

Usage
=
The program takes parameters three parameters such as input file, k, and number of most frequent k-mers to be listed via command line arguments.

* ./kmerCounter "filename" "kmersize" "topcount" 
* Example: ./kmerCounter test.fastq 30 25)

* filename : fastq file name
* kmersize : integer value between 1 and 256
* topcount : integer value between 1 and 10000

* default signature length : 6 (signatureLen = k/2, if (1 < k < 6))

Author : Batuhan Baskurt
Date   : 2017-09-01
Version: 0.21 
