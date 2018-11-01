#include "kmerResultCollector.h"
#include "kmerBinController.h"

kmerCounting::kmerResultCollector* kmerCounting::kmerResultCollector::s_instance = nullptr;
std::mutex kmerCounting::kmerResultCollector::m_instanceMutex;

kmerCounting::kmerResultCollector::kmerResultCollector() :
m_minTop(0),
m_totalNumberOfKmers(0)
{

}


kmerCounting::kmerResultCollector::~kmerResultCollector()
{
}

kmerCounting::kmerResultCollector* kmerCounting::kmerResultCollector::getInstance()
{
	std::lock_guard<std::mutex> lock(m_instanceMutex);

	if (s_instance == NULL)
		s_instance = new kmerResultCollector();
	return s_instance;
}

bool kmerCounting::kmerResultCollector::init(const kmerCounting::kmerParams& params)
{
	m_params = params;

	for (int i = 0; i < m_params.topkMerCount; i++)
	{
		m_topKmers.push_back(std::make_pair("", 0));
	}

	return true;
}

void kmerCounting::kmerResultCollector::notifyResult(const std::vector<std::pair<std::string, int64> > results, int64 totalNumberOfKmersCounted)
{
	std::lock_guard<std::mutex> lock(m_instanceMutex);
	m_totalNumberOfKmers += totalNumberOfKmersCounted;
	std::vector<std::pair<std::string, int64> >::const_iterator it;
	for (it = results.begin(); it != results.end(); it++)
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

void kmerCounting::kmerResultCollector::countingCompleted(long long processTime)
{
	std::lock_guard<std::mutex> lock(m_instanceMutex);
	std::vector<std::pair<std::string, int64> >::const_iterator it;

	std::ofstream fmatch("topKmers.txt", std::ios::out);
	std::cout << "Top kmers: " << std::endl;
	int i = 1;
	for (it = m_topKmers.begin(); it != m_topKmers.end(); it++)
	{
		std::cout << i << ") " << (*it).first << " : " << (*it).second << std::endl;
		fmatch << i++ << ") " << (*it).first << " : " << (*it).second << std::endl;
	}

	std::cout << "Total number of kmers: " << m_totalNumberOfKmers << std::endl;
	fmatch << "Total number of kmers: " << m_totalNumberOfKmers << std::endl;
	fmatch << "Process time: " << processTime << " ms" << std::endl;
	if (m_totalNumberOfKmers != kmerCounting::kmerBinController::getInstance()->getTotalNumberOfKmers())
	{
		std::cout << "Error: Number of extracted kmers and counted kmers are different!" << std::endl;
	}

	fmatch.close();
}

