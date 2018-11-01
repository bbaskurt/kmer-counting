#include "kmerSignature.h"
#include "kmerUtils.h"

uint32 kmerCounting::kmerSignature::norm1[];
uint32 kmerCounting::kmerSignature::norm2[];
uint32 kmerCounting::kmerSignature::norm3[];
uint32 kmerCounting::kmerSignature::norm4[];
uint32 kmerCounting::kmerSignature::norm5[];
uint32 kmerCounting::kmerSignature::norm6[];
uint32 kmerCounting::kmerSignature::norm7[];
uint32 kmerCounting::kmerSignature::norm8[];

kmerCounting::kmerSignature::kmerSignature(const int& signatureLen)
{
	init_norm(norm1, 1);
	init_norm(norm2, 2);
	init_norm(norm3, 3);
	init_norm(norm4, 4);
	init_norm(norm5, 5);
	init_norm(norm6, 6);
	init_norm(norm7, 7);
	init_norm(norm8, 8);

	switch (signatureLen)
	{
	case 1:
		m_norm = norm1;
		break;
	case 2:
		m_norm = norm2;
		break;
	case 3:
		m_norm = norm3;
		break;
	case 4:
		m_norm = norm4;
		break;
	case 5:
		m_norm = norm5;
		break;
	case 6:
		m_norm = norm6;
		break;
	case 7:
		m_norm = norm7;
		break;
	case 8:
		m_norm = norm8;
		break;
	default:
		break;
	}
	m_len = signatureLen;
	m_mask = (1 << m_len * 2) - 1;
	m_str = 0;
}


kmerCounting::kmerSignature::~kmerSignature()
{
}

void kmerCounting::kmerSignature::init_norm(uint32* m_norm, uint32 m_len)
{
	uint32 special = 1 << m_len * 2;
	for (uint32 i = 0; i < special; ++i)
	{
		uint32 rev = kmerCounting::kmerUtils::getRev(i, m_len);
		uint32 str_val = kmerCounting::kmerUtils::isAllowed(i, m_len) ? i : special;
		uint32 rev_val = kmerCounting::kmerUtils::isAllowed(rev, m_len) ? rev : special;
		m_norm[i] = MIN(str_val, rev_val);
	}
}

void kmerCounting::kmerSignature::insert(char* seq)
{
	switch (m_len)
	{
	case 5:
		m_str = (seq[0] << 8) + (seq[1] << 6) + (seq[2] << 4) + (seq[3] << 2) + (seq[4]);
		break;
	case 6:
		m_str = (seq[0] << 10) + (seq[1] << 8) + (seq[2] << 6) + (seq[3] << 4) + (seq[4] << 2) + (seq[5]);
		break;
	case 7:
		m_str = (seq[0] << 12) + (seq[1] << 10) + (seq[2] << 8) + (seq[3] << 6) + (seq[4] << 4) + (seq[5] << 2) + (seq[6]);
		break;
	case 8:
		m_str = (seq[0] << 14) + (seq[1] << 12) + (seq[2] << 10) + (seq[3] << 8) + (seq[4] << 6) + (seq[5] << 4) + (seq[6] << 2) + (seq[7]);
		break;
	default:
		break;
	}

	m_currentVal = m_norm[m_str];
}

void kmerCounting::kmerSignature::insert(uchar symb)
{
	m_str <<= 2;
	m_str += symb;
	m_str &= m_mask;

	m_currentVal = m_norm[m_str];
}

uint32 kmerCounting::kmerSignature::get() const
{
	return m_currentVal;
}

bool kmerCounting::kmerSignature::operator==(const kmerCounting::kmerSignature& x)
{
	return m_currentVal == x.m_currentVal;
}

bool kmerCounting::kmerSignature::operator<(const kmerCounting::kmerSignature& x)
{
	return m_currentVal < x.m_currentVal;
}

void kmerCounting::kmerSignature::clear()
{
	m_str = 0;
}

bool kmerCounting::kmerSignature::operator<=(const kmerCounting::kmerSignature& x)
{
	return m_currentVal <= x.m_currentVal;
}

void kmerCounting::kmerSignature::set(const kmerCounting::kmerSignature& x)
{
	m_str = x.m_str;
	m_currentVal = x.m_currentVal;
}
