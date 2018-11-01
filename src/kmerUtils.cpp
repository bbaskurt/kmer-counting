#include "kmerUtils.h"

kmerCounting::kmerUtils::kmerUtils()
{
}


kmerCounting::kmerUtils::~kmerUtils()
{
}

bool kmerCounting::kmerUtils::isAllowed(uint32 mmer, uint32 len)
{
	if ((mmer & 0x3f) == 0x3f)            // TTT suffix
		return false;
	if ((mmer & 0x3f) == 0x3b)            // TGT suffix
		return false;
	if ((mmer & 0x3c) == 0x3c)            // TG* suffix
		return false;

	for (uint32 j = 0; j < len - 3; ++j)
	{
		if ((mmer & 0xf) == 0)                // AA inside
			return false;
		else
			mmer >>= 2;
	}

	if (mmer == 0)            // AAA prefix
		return false;
	if (mmer == 0x04)        // ACA prefix
		return false;
	if ((mmer & 0xf) == 0)    // *AA prefix
		return false;

	return true;
}

uint32 kmerCounting::kmerUtils::getRev(uint32 mmer, uint32 len)
{
	uint32 rev = 0;
	uint32 shift = len * 2 - 2;
	for (uint32 i = 0; i < len; ++i)
	{
		rev += (3 - (mmer & 3)) << shift;
		mmer >>= 2;
		shift -= 2;
	}
	return rev;
}

char kmerCounting::kmerUtils::getSymbol(char value)
{
	//uint32 x = (m_data[p >> 5] >> (2 * (p & 31))) & 0x03;

	switch (value)
	{
	case 0: return 'A';
	case 1: return 'C';
	case 2: return 'G';
	default: return 'T';
	}
}
