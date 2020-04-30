#include "myutils.h"
#include "linereader.h"

LineReader::LineReader()
	{
	Clear();
	}

LineReader::~LineReader()
	{
	Close();
	}

void LineReader::Open(const string &FileName)
	{
	m_FileName = FileName;
	m_gz = EndsWith(FileName, ".gz");
	if (m_gz)
		m_f = OpenGzipFile(FileName);
	else
		m_f = OpenStdioFile(FileName);
	m_Buffer = myalloc(char, LR_BUFF);
	m_BufferBytes = 0;
	m_BufferOffset = 0;
	m_LineNr = 0;
	m_EOF = false;
	}

void LineReader::Clear()
	{
	m_f = 0;
	m_Buffer = 0;
	m_BufferOffset = 0;
	m_BufferBytes = 0;
	m_LineNr = 0;
	m_gz = false;
	m_EOF = true;
	}

void LineReader::Close()
	{
	if (m_f == 0)
		return;

	if (m_gz)
		CloseGzipFile(m_f);
	else
		CloseStdioFile(m_f);
	m_f = 0;
	myfree(m_Buffer);
	Clear();
	}

bool LineReader::ReadLine(t_LineBuff &Line)
	{
	if (m_EOF)
		return false;

	Line.Alloc(32*1024);
	char *LineData = Line.Data;
	unsigned Length = 0;
	for (;;)
		{
		if (m_BufferOffset >= m_BufferBytes)
			{
			FillBuff();
			if (m_EOF)
				{
				Line.Size = Length;
				if (Length == 0)
					return false;

				LineData[Length] = 0;
				++m_LineNr;
				return true;
				}
			}

		char c = m_Buffer[m_BufferOffset++];
		if (c == '\r')
			continue;
		if (c == '\n')
			{
			LineData[Length] = 0;
			Line.Size = Length;
			++m_LineNr;
			return true;
			}
		if (Length == Line.MaxSize)
			{
			StartTimer(LR_Realloc);
			Line.Size = Length;
			Line.Alloc(Line.MaxSize + 32*1024);
			LineData = Line.Data;
			EndTimer(LR_Realloc);
			}
		LineData[Length++] = c;
		}
	}

void LineReader::FillBuff()
	{
	if (m_EOF)
		return;
	uint32 BytesToRead = LR_BUFF;
	m_BufferOffset = 0;
	if (m_gz)
		m_BufferBytes = ReadGzipFile(m_f, m_Buffer, LR_BUFF);
	else
		m_BufferBytes = ReadStdioFile_NoFail(m_f, m_Buffer, LR_BUFF);
	if (m_BufferBytes == 0)
		m_EOF = true;
	}

void cmd_gunzip()
	{
	const string &InputFileName = opt(gunzip);
	const string &OutputFileName = opt(output);

	FILE *fIn = OpenGzipFile(InputFileName);
	FILE *fOut = CreateStdioFile(OutputFileName);

	const unsigned M = 1024*1024;
	byte *Buffer = myalloc(byte, M);

	for (;;)
		{
		uint32 n = ReadGzipFile(fIn, Buffer, M);
		if (n == 0)
			break;
		WriteStdioFile(fOut, Buffer, n);
		}

	CloseGzipFile(fIn);
	CloseStdioFile(fOut);
	}
