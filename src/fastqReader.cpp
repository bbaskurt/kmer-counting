#include "fastqReader.h"

kmerCounting::fastqReader::fastqReader(bool memoryMode) :
m_memoryMode(memoryMode),
m_filePos(0),
m_fileSize(-1),
m_part(nullptr),
m_partPos(-1),
m_partSize(KMER_FILE_READER_PART_SIZE)
{
	m_part = new uchar[KMER_FILE_READER_PART_SIZE];
	memset(m_part, 0, KMER_FILE_READER_PART_SIZE);
	// Prepare encoding of symbols
	for (int i = 0; i < 256; ++i)
		m_codes[i] = -1;
	m_codes['A'] = m_codes['a'] = 0;
	m_codes['C'] = m_codes['c'] = 1;
	m_codes['G'] = m_codes['g'] = 2;
	m_codes['T'] = m_codes['t'] = 3;
}


kmerCounting::fastqReader::~fastqReader()
{
	if (m_part != nullptr)
	{
		delete[] m_part;
		m_part = nullptr;
	}
}

void kmerCounting::fastqReader::open(const std::string& f_name)
{
	if (m_memoryMode)
	{
		// NOT IMPLEMENTED
	}
	else
	{
		m_file = fopen(f_name.c_str(), "rb"/*"wb+"*/);

		if (!m_file)
		{
			std::cout << "Error: Cannot open temporary file " << f_name << "\n";
			exit(1);
		}
		setbuf(m_file, nullptr);
	}
}

int kmerCounting::fastqReader::close()
{
	if (m_memoryMode)
	{
		// NOT IMPLEMENTED
		return 0;
	}
	else
	{
		return fclose(m_file);
	}
}

size_t kmerCounting::fastqReader::read(uchar * ptr, size_t size, size_t count)
{
	if (m_memoryMode)
	{
		uint64 pos = 0;
		// NOT IMPLEMENTED
		return pos;
	}
	else
	{
		return fread(ptr, size, count, m_file);
	}
}

size_t kmerCounting::fastqReader::write(const uchar * ptr, size_t size, size_t count)
{
	if (m_memoryMode)
	{
		// NOT IMPLEMENTED
		return size * count;
	}
	else
	{
		size_t readData = fwrite(ptr, size, count, m_file);
		m_filePos += readData;
		return readData;
	}
}

bool kmerCounting::fastqReader::getSequence(char *seq, uint32 &seqSize)
{
	uchar c = 0;
	uint32 pos = 0;
	//std::cout << "Start" << std::endl;

	if (!isPartAvailable())
	{
		return false;
	}


	c = m_part[m_partPos++];

	if (c != '@')
	{
		//std::cout << "Error on getting sequence" << std::endl;
		return false;
	}

	for (; m_partPos < m_partSize;)
	{
		c = m_part[m_partPos++];
		if (c < 32)					// newliners
			break;
	}
	if (!isPartAvailable())
	{
		return false;
	}
		

	c = m_part[m_partPos++];
	if (c >= 32)
	{
		m_partPos--;
	}
	else if (m_partPos >= m_partSize)
	{
		if (!isPartAvailable())
		{
			return false;
		}
	}

	// Sequence
	for (; m_partPos < m_partSize;)
	{
		c = m_part[m_partPos++];
		if (c < 32)					// newliners
			break;
		seq[pos++] = m_codes[c];
	}

	if (!isPartAvailable())
	{
		return false;
	}

	c = m_part[m_partPos++];
	if (c >= 32)
		m_partPos--;
	else if (m_partPos >= m_partSize)
	{
		if (!isPartAvailable())
		{
			return false;
		}
	}

	// Plus
	c = m_part[m_partPos++];
	if (m_partPos >= m_partSize)
	{
		if (!isPartAvailable())
		{
			return false;
		}
	}
	if (c != '+')
		return false;
	for (; m_partPos < m_partSize;)
	{
		c = m_part[m_partPos++];
		if (c < 32)					// newliners
			break;
	}
	if (!isPartAvailable())
	{
		return false;
	}

	c = m_part[m_partPos++];
	if (c >= 32)
		m_partPos--;
	else if (m_partPos >= m_partSize)
	{
		if (!isPartAvailable())
		{
			return false;
		}			
	}

	// Quality
	m_partPos += pos;
	if (!isPartAvailable())
	{
		return false;
	}
	c = m_part[m_partPos++];
	seqSize = pos;

	if (m_partPos >= m_partSize)
		return true;

	if (m_part[m_partPos++] >= 32)
		m_partPos--;
	else if (m_partPos >= m_partSize)
	{
		return true;
	}

	//std::cout << "End of get sequence: " << m_partPos << "-" << (c == '\n' || c == '\r') << std::endl;
	return (c == '\n' || c == '\r');
}

bool kmerCounting::fastqReader::isPartAvailable()
{
	if (m_partPos >= m_partSize || m_partPos == -1)
	{
		m_partSize = read(m_part, 1, KMER_FILE_READER_PART_SIZE);
		if (m_partSize > 0)
		{
			m_filePos += m_partSize;
			m_partPos = 0;
		}
		else
		{
			if (feof(m_file) == 0)		// if it is not eof
			{
				std::cout << "Error on getting sequence" << std::endl;
				return false;
			}
		}
	}
	return true;
}