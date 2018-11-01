#include "kmerExtractor.h"
#include "kmerUtils.h"
#include "kmerSignature.h"
#include "fastqReader.h"
#include "kmerBinController.h"

kmerCounting::kmerExtractor::kmerExtractor(const kmerCounting::kmerParams& params) :
m_plusXRecsCount(0)
{
	m_params = params;
	m_fastqReader = new fastqReader(false);
	m_fastqReader->open(params.filename);

	// total number of possible signatures for given signature length
	m_mapSize = (1 << 2 * m_params.signatureLen) + 1;
	m_signatureMap = new int32[m_mapSize];
	std::fill_n(m_signatureMap, m_mapSize, -1);

	// special signature to be used for kmers that does not have any valid signature
	uint32 special = 1 << m_params.signatureLen * 2;
	for (int32 i = 0; i < m_mapSize; ++i)
	{
		uint32 rev = kmerCounting::kmerUtils::getRev(i, m_params.signatureLen);
		uint32 str_val = kmerCounting::kmerUtils::isAllowed(i, m_params.signatureLen) ? i : special;
		uint32 rev_val = kmerCounting::kmerUtils::isAllowed(rev, m_params.signatureLen) ? rev : special;
		uint32 sigVal = MIN(str_val, rev_val);
		m_signatureMap[i] = sigVal;
	}
}

kmerCounting::kmerExtractor::~kmerExtractor()
{
	if (m_fastqReader)
	{
		//m_fastqReader->close();
		delete m_fastqReader;
		m_fastqReader = nullptr;
	}
}

void kmerCounting::kmerExtractor::putExtendedKmer(uint32 binNo, char* seq, uint32 n)
{
	//// TEST
	//for (int i = 0; i < n; i++)
	//{
	//	std::cout << (int)seq[i] << "-";
	//}
	//std::cout << std::endl;

	// allocate buffer if current bin is being operated first time
	std::unordered_map<int32, std::pair<uchar*, int32> >::iterator it;
	it = m_bins.find(binNo);
	if (it == m_bins.end())
	{
		// delete if file exists
		std::string tmpFileName = getFileName(binNo, -1);
		std::string fullFilename = tmpFileName + ".bin";
		std::remove(fullFilename.c_str());

		uchar* buff = new uchar[KMER_BIN_BUFFER_SIZE];
		m_bins[binNo] = std::pair<uchar*, int32>(buff, 0);
	}

	// write data to file if buffer size exceed to upper limit.
	uint32 bytes = 1 + (n + 3) / 4;
	if (m_bins[binNo].second + n > KMER_BIN_BUFFER_SIZE)
	{
		//write current buffer to file
		std::string tmpFileName = getFileName(binNo, n);
		std::string fullFilename = tmpFileName + ".bin";
		FILE* file = fopen(fullFilename.c_str(), "ab"/*"wb+"*/);
		size_t dataWritten = fwrite(m_bins[binNo].first, sizeof(char), m_bins[binNo].second, file);
		fclose(file);

		// kmer and super kmer count of the file is updated with local counts also. 
		// update operation must be done using file name taken from getFileName, because 
		// it might be supplementary file in case of main file of signature is full.
		kmerCounting::kmerBinController::getInstance()->addBinSize(tmpFileName, (int32)dataWritten);
		kmerCounting::kmerBinController::getInstance()->addKmerCountToBin(tmpFileName, m_kmerCount[binNo]);
		kmerCounting::kmerBinController::getInstance()->addSuperKmerCountToBin(tmpFileName, m_superKmersCount[binNo]);
		m_kmerCount[binNo] = 0;
		m_superKmersCount[binNo] = 0;
		// set buffer pointer to initial position
		m_bins[binNo].second = 0;
	}

	// TODO: bitwise operations improve the performance
	//// bitwise operations to store each symbol in 2 bits.
	//m_bins[binNo].first[m_bins[binNo].second++] = n - m_params.kMerSize;
	////buffer[buffer_pos++] = n - kmer_len;
	//// write 4 bytes in each iteration by bitwise operation, then remaining part of mod4
	//for (uint32 i = 0, j = 0; i < n / 4; ++i, j += 4)
	//	m_bins[binNo].first[m_bins[binNo].second++] = (seq[j] << 6) + (seq[j + 1] << 4) + (seq[j + 2] << 2) + seq[j + 3];
	//switch (n % 4)
	//{
	//case 1:
	//	m_bins[binNo].first[m_bins[binNo].second++] = (seq[n - 1] << 6);
	//	break;
	//case 2:
	//	m_bins[binNo].first[m_bins[binNo].second++] = (seq[n - 2] << 6) + (seq[n - 1] << 4);
	//	break;
	//case 3:
	//	m_bins[binNo].first[m_bins[binNo].second++] = (seq[n - 3] << 6) + (seq[n - 2] << 4) + (seq[n - 1] << 2);
	//	break;
	//}

	// write number of symbols
	m_bins[binNo].first[m_bins[binNo].second++] = n;
	// write each symbol in 1-byte
	for (uint32 i = 0; i < n; i++)
	{
		// TEST
		//std::cout << (char)kmerCounting::kmerUtils::getSymbol(seq[i]);	// << std::endl;
		m_bins[binNo].first[m_bins[binNo].second++] = seq[i];
	}
	//std::cout << std::endl;

	increaseSuperKmerCount(binNo, 1);
	increaseKmerCount(binNo, n - m_params.kMerSize + 1);
}


void kmerCounting::kmerExtractor::extractSuperkmers()
{
	char* seq = new char[KMER_MAX_SEQUENCE_SIZE];
	uint32 seq_size;

	uint32 signature_start_pos;
	kmerSignature current_signature(m_params.signatureLen), end_mmer(m_params.signatureLen);
	uint32 bin_no;

	uint32 i;
	int len;		//length of extended kmer

	while (m_fastqReader->getSequence(seq, seq_size))
	{
		i = 0;
		len = 0;
		while (i + m_params.kMerSize - 1 < seq_size)
		{
			bool contains_N = false;
			//building first kmerSignature after 'N' or at the read begining
			for (int j = 0; j < m_params.signatureLen; ++j, ++i)
				if (seq[i] < 0)//'N'
				{
					contains_N = true;
					break;
				}
			//kmerSignature must be shorter than k-mer so if kmerSignature contains 'N', k-mer will contains it also
			if (contains_N)
			{
				++i;
				continue;
			}
			len = m_params.signatureLen;
			signature_start_pos = i - m_params.signatureLen;
			current_signature.insert(seq + signature_start_pos);
			end_mmer.set(current_signature);
			for (; i < seq_size; ++i)
			{
				if (seq[i] < 0)//'N'
				{
					if (len >= m_params.kMerSize)
					{
						bin_no = m_signatureMap[current_signature.get()];
						putExtendedKmer(bin_no, seq + i - len, len);
					}
					len = 0;
					++i;
					break;
				}
				end_mmer.insert(seq[i]);
				if (end_mmer < current_signature)		//kmerSignature at the end of current k-mer is lower than current
				{
					if (len >= m_params.kMerSize)
					{
						bin_no = m_signatureMap[current_signature.get()];
						putExtendedKmer(bin_no, seq + i - len, len);
						len = m_params.kMerSize - 1;
					}
					current_signature.set(end_mmer);
					signature_start_pos = i - m_params.signatureLen + 1;
				}
				else if (end_mmer == current_signature)
				{
					current_signature.set(end_mmer);
					signature_start_pos = i - m_params.signatureLen + 1;
				}
				else if (signature_start_pos + m_params.kMerSize - 1 < i)		//need to find new kmerSignature
				{
					bin_no = m_signatureMap[current_signature.get()];
					putExtendedKmer(bin_no, seq + i - len, len);
					len = m_params.kMerSize - 1;
					//looking for new kmerSignature
					++signature_start_pos;
					//building first kmerSignature in current k-mer
					end_mmer.insert(seq + signature_start_pos);
					current_signature.set(end_mmer);
					for (uint32 j = signature_start_pos + m_params.signatureLen; j <= i; ++j)
					{
						end_mmer.insert(seq[j]);
						if (end_mmer <= current_signature)
						{
							current_signature.set(end_mmer);
							signature_start_pos = j - m_params.signatureLen + 1;
						}
					}
				}
				++len;
				if (len == m_params.kMerSize + 255)			//one byte is used to store counter of additional symbols in extended k-mer
				{
					bin_no = m_signatureMap[current_signature.get()];
					putExtendedKmer(bin_no, seq + i + 1 - len, len);
					i -= m_params.kMerSize - 2;
					len = 0;
					break;
				}

			}
		}
		if (len >= m_params.kMerSize)		//last one in read
		{
			bin_no = m_signatureMap[current_signature.get()];
			putExtendedKmer(bin_no, seq + i - len, len);
		}
	}
	delete[] seq;

	//write remaining buffer to file
	std::unordered_map<int32, std::pair<uchar*, int32> >::iterator it;
	for (it = m_bins.begin(); it != m_bins.end(); it++)
	{
		// TODO: merge buffers according to containing number of kmers in order to reduce I/O cost
		// Be careful about not merging supplementary files. All instances of a kmer must be in same bin + supplementary file.

		std::string tmpFileName = getFileName((*it).first, (*it).second.second);
		std::string fullFilename = tmpFileName + ".bin";
		FILE* file = fopen(fullFilename.c_str(), "ab"/*"wb+"*/);
		size_t dataWritten = fwrite((*it).second.first, sizeof(char), (*it).second.second, file);
		fclose(file);
		kmerCounting::kmerBinController::getInstance()->addBinSize(tmpFileName, (int32)dataWritten);
		kmerCounting::kmerBinController::getInstance()->addKmerCountToBin(tmpFileName, m_kmerCount[(*it).first]);
		kmerCounting::kmerBinController::getInstance()->addSuperKmerCountToBin(tmpFileName, m_superKmersCount[(*it).first]);
		m_kmerCount[(*it).first] = 0;
		m_superKmersCount[(*it).first] = 0;
		if ((*it).second.first != nullptr)
		{
			delete[](*it).second.first;
			(*it).second.first = nullptr;
		}
	}

	m_fastqReader->close();
}

std::string kmerCounting::kmerExtractor::getFileName(uint32 binNo, int32 bytes)
{
	std::string tmpFileName = std::to_string(binNo);

	while (tmpFileName.length() < 5)
		tmpFileName = std::string("0") + tmpFileName;

	if (bytes != -1)
	{
		// Calculate estimated file size after current write. If it is bigger than maximum allowed size, 
		// create new file for the signature and register it to bin controller.
		int64 fileSize = kmerCounting::kmerBinController::getInstance()->getBinSize(tmpFileName);
		if ((fileSize + m_bins[binNo].second + bytes) > KMER_MAX_TEMPORARY_FILE_SIZE)
		{
			std::vector<std::string> fileList = kmerCounting::kmerBinController::getInstance()->getFileListOfBin(std::to_string(binNo));
			bool createNewTmpFile = false;
			if (fileList.empty())
			{
				createNewTmpFile = true;
			}
			else
			{
				fileSize = kmerCounting::kmerBinController::getInstance()->getBinSize(fileList.at(fileList.size() - 1));
				if ((fileSize + m_bins[binNo].second + bytes) > KMER_MAX_TEMPORARY_FILE_SIZE)
				{
					createNewTmpFile = true;
				}
				else
				{
					tmpFileName = fileList.at(fileList.size() - 1);
				}
			}
			if (createNewTmpFile)
			{
				tmpFileName = tmpFileName + "_" + std::to_string(fileList.size() + 1);
				kmerCounting::kmerBinController::getInstance()->addFileToBin(std::to_string(binNo), tmpFileName);
				// delete if any old file exists with same name. (from previous execution)
				std::remove(tmpFileName.c_str());
			}
		}

	}

	return tmpFileName;
}

uint32 kmerCounting::kmerExtractor::increaseSuperKmerCount(uint32 binNo, uint32 increment)
{
	std::unordered_map<int32, uint32>::iterator itCount;
	itCount = m_superKmersCount.find(binNo);
	if (itCount != m_superKmersCount.end())
	{
		m_superKmersCount[binNo] += increment;
	}
	else
	{
		m_superKmersCount[binNo] = increment;
	}

	return m_superKmersCount[binNo];
}


uint32 kmerCounting::kmerExtractor::increaseKmerCount(uint32 binNo, uint32 increment)
{
	std::unordered_map<int32, uint32>::iterator itCount;
	itCount = m_kmerCount.find(binNo);
	if (itCount != m_kmerCount.end())
	{
		m_kmerCount[binNo] += increment;
	}
	else
	{
		m_kmerCount[binNo] = increment;
	}

	return m_kmerCount[binNo];
}