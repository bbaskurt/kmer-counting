#include "kmc.h"
#include "fastqReader.h"
#include "kmerSignature.h"
#include "kmerUtils.h"
#include "kmerExtractor.h"
#include "kmerCounter.h"
#include "kmerBinController.h"
#include "kmerResultCollector.h"

kmerCounting::kmc::kmc(const kmerCounting::kmerParams& params)
{
	auto startTime = std::chrono::system_clock::now();
	m_params = params;

	// TODO: one thread reads input file and pushs to a queue, other threads process them parallel.
	std::cout << "kmer extraction in progress ..." << '\n';
	m_kmerExtractor = new kmerCounting::kmerExtractor(params);

	m_kmerExtractor->extractSuperkmers();

	delete m_kmerExtractor;
	m_kmerExtractor = nullptr;

	auto endTime = std::chrono::system_clock::now();
	auto dur = endTime - startTime;
	auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(dur);
	std::cout << "kmer extraction completed in: " << ms.count() << " ms" << '\n';
	std::cout << "kmer counting in progress ..." << '\n';

	// kmer counting must wait until all kmers are extracted.

	kmerCounting::kmerResultCollector::getInstance()->init(m_params);

	std::vector<int32> binIdList;
	kmerCounting::kmerBinController::getInstance()->getAllBinIds(binIdList);

	std::vector<kmerCounter*> counters;
	std::vector<std::thread> counterThreads;
	int binCountForThread = (int)binIdList.size() / m_params.counterThreadCount;
	if (binIdList.size() % m_params.counterThreadCount != 0)
	{
		binCountForThread++;
	}
	int binOffset = 0;
	for (int i = 0; i < m_params.counterThreadCount; ++i)
	{
		kmerCounter *counter = new kmerCounter(params);
		counters.push_back(counter);
		std::vector<int32> binIds;
		if (binOffset + binCountForThread >= (int)binIdList.size())
			binCountForThread = (int)binIdList.size() - binOffset;
		for (int b = binOffset; b < binOffset+binCountForThread; b++)
		{
			binIds.push_back(binIdList[b]);
		}
		binOffset += binCountForThread;
		counterThreads.push_back(std::thread(&kmerCounter::processBins, counter, binIds));
		//// TODO: balance work load of threads
		//// TEST
		//int totalSizeTemp = 0;
		//for (int tt = 0; tt < binIds.size(); tt++)
		//{
		//	std::string tmpFileName = std::to_string(binIds[tt]);

		//	while (tmpFileName.length() < 5)
		//		tmpFileName = std::string("0") + tmpFileName;
		//	totalSizeTemp += kmerCounting::kmerBinController::getInstance()->getBinSize(tmpFileName);
		//}
		//std::cout << "Total file size for thread (" << i << "): " << totalSizeTemp << std::endl;
		//// END TEST
	}

	for (auto p = counterThreads.begin(); p != counterThreads.end(); ++p)
		p->join();

	for (auto p = counters.begin(); p != counters.end(); ++p)
		delete *p;

	//std::cout << "End of kmc" << std::endl;
	endTime = std::chrono::system_clock::now();
	dur = endTime - startTime;
	ms = std::chrono::duration_cast<std::chrono::milliseconds>(dur);
	std::cout << "Process time: " << ms.count() << " ms" << '\n';

	kmerCounting::kmerResultCollector::getInstance()->countingCompleted(ms.count());

	kmerCounting::kmerBinController::getInstance()->removeTempFiles();
}

kmerCounting::kmc::~kmc()
{
	if (m_kmerExtractor)
	{
		delete m_kmerExtractor;
		m_kmerExtractor = nullptr;
	}
}

