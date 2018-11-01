#include "kmerCounter.h"
#include "kmerBaseTypes.h"
#include "kmerBinController.h"
#include <stddef.h>
#include "kmerUtils.h"
#include "kmerResultCollector.h"

kmerCounting::kmerCounter::kmerCounter(const kmerCounting::kmerParams& params):
m_minTop(0)
{
	m_params = params;
	m_maxKmerSizeForBin = kmerCounting::kmerBinController::getInstance()->getMaxKmerCountForBin();

	for (int i = 0; i < m_params.topkMerCount; i++)
	{
		m_topKmers.push_back(std::make_pair("", 0));
	}
}


kmerCounting::kmerCounter::~kmerCounter()
{

}

void kmerCounting::kmerCounter::updateTopKmers()
{
	// update top kmers
	std::unordered_map<std::string, int64>::iterator it;
	for (it = m_kmerCounts.begin(); it != m_kmerCounts.end(); it++)
	{
		if ((*it).second > m_minTop)
		{
			for (int i = 0; i < m_params.topkMerCount; i++)
			{
				if (m_topKmers[i].second < (*it).second)
				{
					m_topKmers.insert(m_topKmers.begin() + i, std::make_pair((*it).first, (*it).second));
					m_topKmers.erase(m_topKmers.begin() + m_topKmers.size() - 1);
					m_minTop = m_topKmers[m_topKmers.size() - 1].second;
					break;
				}
			}
		}
	}
}

void kmerCounting::kmerCounter::processBins(const std::vector<int32> binIds)
{
	int64 totalCount = 0;
	for (auto binId : binIds)
	{
		//binsChecker.push_back(binId);
		m_kmerCounts.clear();
		processBin(binId);
		updateTopKmers();

		std::unordered_map<std::string, int64>::iterator it;
		for (it = m_kmerCounts.begin(); it != m_kmerCounts.end(); it++)
		{
			totalCount += (*it).second;
		}
	}

	kmerCounting::kmerResultCollector::getInstance()->notifyResult(m_topKmers, totalCount);
}

void kmerCounting::kmerCounter::processBin(const int32& binNo)
{
	// read bin
	std::vector<std::string> fileList = kmerCounting::kmerBinController::getInstance()->getFileListOfBin(std::to_string(binNo));
	// add original file to list
	std::string tmpFileName = std::to_string(binNo);
	while (tmpFileName.length() < 5)
		tmpFileName = std::string("0") + tmpFileName;
	fileList.push_back(tmpFileName);

	for (auto filename : fileList)
	{
		uchar* bin = new uchar[kmerCounting::kmerBinController::getInstance()->getBinSize(filename)];
		memset(bin, 0, kmerCounting::kmerBinController::getInstance()->getBinSize(filename));
		std::string fullFileName = filename + ".bin";
		FILE* file = fopen(fullFileName.c_str(), "rb"/*"wb+"*/);
		uint64 readed = fread(bin, 1/*sizeof(uint32)*/, kmerCounting::kmerBinController::getInstance()->getBinSize(filename), file);
		//uint32 kmerCount = kmerCounting::kmerBinController::getInstance()->getKmerCount(filename);
		fclose(file);

		uint64 i = 0;
		while (i < readed)
		{
			int superKmerLength = bin[i];
			for (int j = 1; j < superKmerLength - m_params.kMerSize + 2; j++)
			{
				std::string kmerRead = "";
				for (int k = 0; k < m_params.kMerSize; k++)
				{
					kmerRead += kmerCounting::kmerUtils::getSymbol(bin[i+j+k]);
				}			
				countKmer(kmerRead);
			}
			i += superKmerLength+1;
		}

		delete[] bin;
	}
}

void kmerCounting::kmerCounter::countKmer(const std::string& kmer)
{
	std::unordered_map<std::string, int64>::iterator it;
	it = m_kmerCounts.find(kmer);
	if (it != m_kmerCounts.end())
	{
		m_kmerCounts[kmer]++;
	}
	else
	{
		m_kmerCounts[kmer] = 1;
	}
}