#include "kmerBinController.h"

kmerCounting::kmerBinController* kmerCounting::kmerBinController::s_instance = nullptr;
std::mutex kmerCounting::kmerBinController::m_instanceMutex;

kmerCounting::kmerBinController::kmerBinController():
m_maxKmerCountForBin(0)
{
}


kmerCounting::kmerBinController::~kmerBinController()
{
}

kmerCounting::kmerBinController* kmerCounting::kmerBinController::getInstance()
{
	std::lock_guard<std::mutex> lock(m_instanceMutex);

	if (s_instance == NULL)
		s_instance = new kmerBinController();
	return s_instance;
}

void kmerCounting::kmerBinController::addBinSize(std::string binId, int64 binSize)
{
	std::lock_guard<std::mutex> lock(m_instanceMutex);
	std::unordered_map<std::string, int64>::iterator it;
	it = m_binSizeMap.find(binId);
	if (it == m_binSizeMap.end())
	{
		m_binSizeMap[binId] = binSize;
	}
	else
	{
		m_binSizeMap[binId] += binSize;
	}
}

int64 kmerCounting::kmerBinController::getBinSize(std::string binId)
{
	std::lock_guard<std::mutex> lock(m_instanceMutex);
	std::unordered_map<std::string, int64>::iterator it;
	it = m_binSizeMap.find(binId);
	if (it == m_binSizeMap.end())
	{
		return 0;
	}
	return m_binSizeMap[binId];
}

int kmerCounting::kmerBinController::addFileToBin(std::string binId, const std::string& filename)
{
	std::lock_guard<std::mutex> lock(m_instanceMutex);
	std::unordered_map<std::string, std::vector<std::string> >::iterator it;
	it = m_binsWithMultipleFiles.find(binId);
	if (it == m_binsWithMultipleFiles.end())
	{
		std::vector<std::string> fileList;
		fileList.push_back(filename);
		m_binsWithMultipleFiles[binId] = fileList;
		return 1;
	}
	else
	{
		(*it).second.push_back(filename);
		return (int)(*it).second.size();
	}
	return -1;
}

std::vector<std::string> kmerCounting::kmerBinController::getFileListOfBin(std::string binId)
{
	std::lock_guard<std::mutex> lock(m_instanceMutex);
	std::unordered_map<std::string, std::vector<std::string> >::iterator it;
	it = m_binsWithMultipleFiles.find(binId);
	if (it != m_binsWithMultipleFiles.end())
	{
		return (*it).second;
	}
	else
	{
		return{};
	}
	return{};
}

void kmerCounting::kmerBinController::getAllBinIds(std::vector<int32>& binIdList)
{
	std::lock_guard<std::mutex> lock(m_instanceMutex);
	std::unordered_map<std::string, int64>::iterator it;
	binIdList.clear();

	for (it = m_binSizeMap.begin(); it != m_binSizeMap.end(); it++)
	{
		std::size_t found = (*it).first.find("_");
		if (found == std::string::npos)
		{
			binIdList.push_back(std::stoi((*it).first));
		}
	}
}

void kmerCounting::kmerBinController::addKmerCountToBin(std::string binId, int64 kmerCount)
{
	std::lock_guard<std::mutex> lock(m_instanceMutex);
	std::unordered_map<std::string, int64>::iterator it;
	it = m_kmerCountForBin.find(binId);
	if (it == m_kmerCountForBin.end())
	{
		m_kmerCountForBin[binId] = kmerCount;
	}
	else
	{
		m_kmerCountForBin[binId] += kmerCount;
	}
	if (m_kmerCountForBin[binId] > m_maxKmerCountForBin)
		m_maxKmerCountForBin = m_kmerCountForBin[binId];
}

int64 kmerCounting::kmerBinController::getKmerCount(std::string binId)
{
	std::lock_guard<std::mutex> lock(m_instanceMutex);
	//std::vector<std::string> fileList = getFileListOfBin(binId);
	//int64 totalCount = 0;
	//fileList.push_back(binId);
	//for (auto file : fileList)
	//{
	//	std::unordered_map<std::string, int64>::iterator it;
	//	it = m_kmerCountForBin.find(file);
	//	if (it != m_kmerCountForBin.end())
	//	{
	//		totalCount += m_kmerCountForBin[file];
	//	}
	//}
	//return totalCount;

	std::unordered_map<std::string, int64>::iterator it;
	it = m_kmerCountForBin.find(binId);
	if (it == m_kmerCountForBin.end())
	{
		return -1;
	}
	return m_kmerCountForBin[binId];
	
}

void kmerCounting::kmerBinController::addSuperKmerCountToBin(std::string binId, int64 superKmerCount)
{
	std::lock_guard<std::mutex> lock(m_instanceMutex);
	std::unordered_map<std::string, int64>::iterator it;
	it = m_superkmerCountForBin.find(binId);
	if (it == m_superkmerCountForBin.end())
	{
		m_superkmerCountForBin[binId] = superKmerCount;
	}
	else
	{
		m_superkmerCountForBin[binId] += superKmerCount;
	}
}

int64 kmerCounting::kmerBinController::getSuperKmerCount(std::string binId)
{
	std::lock_guard<std::mutex> lock(m_instanceMutex);
	std::unordered_map<std::string, int64>::iterator it;
	it = m_superkmerCountForBin.find(binId);
	if (it == m_superkmerCountForBin.end())
	{
		return -1;
	}
	return m_superkmerCountForBin[binId];

}

int64 kmerCounting::kmerBinController::getTotalNumberOfKmers()
{
	std::lock_guard<std::mutex> lock(m_instanceMutex);
	int64 totalCount = 0;
	std::unordered_map<std::string, int64>::iterator it;
	for (it = m_kmerCountForBin.begin(); it != m_kmerCountForBin.end(); it++)
	{
		totalCount += (*it).second;
	}
	return totalCount;
}

int64 kmerCounting::kmerBinController::getMaxKmerCountForBin()
{
	std::lock_guard<std::mutex> lock(m_instanceMutex);
	return m_maxKmerCountForBin;
}

void kmerCounting::kmerBinController::removeTempFiles()
{
	std::lock_guard<std::mutex> lock(m_instanceMutex);
	// delete main files
	std::unordered_map<std::string, int64>::iterator it;
	for (it = m_binSizeMap.begin(); it != m_binSizeMap.end(); it++)
	{
		std::string fullName = (*it).first;
		while (fullName.length() < 5)
			fullName = std::string("0") + fullName;
		fullName += ".bin";
		std::remove(fullName.c_str());
	}
	// delete supplementary files
	std::unordered_map<std::string, std::vector<std::string> >::iterator itSupp;
	for (itSupp = m_binsWithMultipleFiles.begin(); itSupp != m_binsWithMultipleFiles.end(); itSupp++)
	{
		for (auto file : (*itSupp).second)
		{
			while (file.length() < 5)
				file = std::string("0") + file;
			file += ".bin";
			std::remove(file.c_str());
		}
	}
	
}