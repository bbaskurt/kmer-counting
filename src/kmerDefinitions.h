#ifndef _KMER_DEFINITIONS_H
#define _KMER_DEFINITIONS_H

#define _CRT_SECURE_NO_WARNINGS
#define _SCL_SECURE_NO_WARNINGS

#include <string>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <memory>
#include <chrono>
#include <algorithm>  
#include <thread> 
#include <time.h>

#define MIN(x,y)	((x) < (y) ? (x) : (y))
#define MAX(x,y)	((x) > (y) ? (x) : (y))
#define NORM(x, lower, upper)	((x) < (lower) ? (lower) : (x) > (upper) ? (upper) : (x))

#define uchar	unsigned char

// partition size for each read of fastq file
#ifndef KMER_FILE_READER_PART_SIZE
#define KMER_FILE_READER_PART_SIZE 104857600	// 100MB
#endif

// super k-mers up to this size are buffered before writing temporary file
#ifndef KMER_BIN_BUFFER_SIZE
#define KMER_BIN_BUFFER_SIZE 256000	// 256K
// SUPPLEMENTARY FILE TEST
//#define KMER_BIN_BUFFER_SIZE 1024	
#endif

// maximum number of characters that can be exist in a line of input file.
#ifndef KMER_MAX_SEQUENCE_SIZE
#define KMER_MAX_SEQUENCE_SIZE 2048	
#endif

// maximum length of k-mer
#ifndef KMER_MAX_K
#define KMER_MAX_K 256	
#endif

// maximum size of temporary file
#ifndef KMER_MAX_TEMPORARY_FILE_SIZE
#define KMER_MAX_TEMPORARY_FILE_SIZE 104857600	// 100 MB	
// SUPPLEMENTARY FILE TEST
//#define KMER_MAX_TEMPORARY_FILE_SIZE 10240
#endif

// maximum number of most frequent kmers to be list.
#ifndef KMER_MAX_TOP_COUNT
#define KMER_MAX_TOP_COUNT 10000
#endif

// maximum number of symbols to be used in coding (A, T, C, G, etc.)
#ifndef MAX_NUMBER_OF_SYMBOLS
#define MAX_NUMBER_OF_SYMBOLS 256	
#endif

// state for finding canonical form of k-mer
enum comparision_state  { kmer_smaller, rev_smaller, equals };

#ifdef WIN32
#define my_fopen	fopen
#define my_fseek	_fseeki64
#define my_ftell	_ftelli64
typedef int int32;
typedef unsigned int uint32;
typedef long long int64;
typedef unsigned long long uint64;

#else
#define my_fopen	fopen
#define my_fseek	fseek
#define my_ftell	ftell
#define _TCHAR	char
#define _tmain	main

typedef int int32;
typedef unsigned int uint32;
typedef long long int64;
typedef unsigned long long uint64;

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#endif

#endif	// _KMER_DEFINITIONS_H
