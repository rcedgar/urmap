#include "myutils.h"
#include "cigar.h"

void PathToCIGAR(const char *Path, string &CIGAR)
	{
	char LastC = *Path;
	unsigned n = 1;
	char Tmp[32];
	for (unsigned i = 1; ; ++i)
		{
		char c = Path[i];
		if (c == 0)
			break;
		if (c == LastC)
			{
			++n;
			continue;
			}
		else
			{
			assert(n > 0);
			if (LastC == 'D')
				LastC = 'I';
			else if (LastC == 'I')
				LastC = 'D';
			sprintf(Tmp, "%u%c", n, LastC);
			CIGAR += string(Tmp);
			LastC = c;
			n = 1;
			}
		}
	if (n > 0)
		{
		if (LastC == 'D')
			LastC = 'I';
		else if (LastC == 'I')
			LastC = 'D';
		sprintf(Tmp, "%u%c", n, LastC);
		CIGAR += string(Tmp);
		}
	}

void CIGARGetOps(const string &CIGAR, vector<char> &Ops,
  vector<unsigned> &Lengths)
	{
	Ops.clear();
	Lengths.clear();
	if (CIGAR.empty())
		return;

	unsigned L = SIZE(CIGAR);
	unsigned n = 0;
	for (unsigned i = 0; i < L; ++i)
		{
		char c = CIGAR[i];
		if (isdigit(c))
			n = n*10 + (c - '0');
		else if (isupper(c) || c == '=')
			{
			if (n == 0)
				Die("Operation '%c' has zero length in CIGAR '%s'", c, CIGAR.c_str());
			Ops.push_back(c);
			Lengths.push_back(n);
			n = 0;
			}
		else
			Die("Invalid char '%c' in CIGAR '%s'", c, CIGAR.c_str());
		}
	if (n > 0)
		Die("Missing operation at end of CIGAR '%s'", CIGAR.c_str());
	}

unsigned CIGARToQL(const string &CIGAR)
	{
	vector<char> Ops;
	vector<unsigned> Lengths;
	CIGARGetOps(CIGAR, Ops, Lengths);
	const unsigned N = SIZE(Ops);
	asserta(SIZE(Lengths) == N);
	unsigned QL = 0;
	for (unsigned i = 0; i < N; ++i)
		{
		char Op = Ops[i];
		switch (Op)
			{
		case 'M':
		case 'I':
			QL += Lengths[i];
			break;

		case 'D':
			break;

		default:
			Die("Unsupported op '%c' in CIGAR '%s'", Op, CIGAR.c_str());
			}
		}
	return QL;
	}

void CIGAROpsToLs(const vector<char> &Ops, const vector<unsigned> &Lengths,
  unsigned &QL, unsigned &TL)
	{
	QL = 0;
	TL = 0;
	const unsigned N = SIZE(Ops);
	asserta(SIZE(Lengths) == N);
	for (unsigned i = 0; i < N; ++i)
		{
		char Op = Ops[i];
		switch (Op)
			{
		case 'M':
			QL += Lengths[i];
			TL += Lengths[i];
			break;

	// CIGAR D&I reverse of my usual convention
		case 'I':
			QL += Lengths[i];
			break;

		case 'D':
			TL += Lengths[i];
			break;

		default:
			Die("Unsupported op '%c' in CIGAR", Op);
			}
		}
	}

void CIGARToLs(const string &CIGAR, unsigned &QL, unsigned &TL)
	{
	vector<char> Ops;
	vector<unsigned> Lengths;
	CIGARGetOps(CIGAR, Ops, Lengths);
	CIGAROpsToLs(Ops, Lengths, QL, TL);
	}

void CIGAROpsFixDanglingMs(vector<char> &Ops, vector<unsigned> &Lengths)
	{
	const unsigned N = SIZE(Ops);
	asserta(SIZE(Lengths) == N);
	if (N < 3)
		return;

// 1M 6I 100M
	if (Ops[0] == 'M' && Lengths[0] <= 2 && Lengths[1] > 4 && Ops[2] == 'M')
		{
		unsigned OldQL;
		unsigned OldTL;
		CIGAROpsToLs(Ops, Lengths, OldQL, OldTL);

		vector<char> NewOps;
		vector<unsigned> NewLengths;
		for (unsigned i = 1; i < N; ++i)
			{
			NewOps.push_back(Ops[i]);
			NewLengths.push_back(Lengths[i]);
			}
		NewLengths[1] += Lengths[0];

		unsigned NewQL;
		unsigned NewTL;
		CIGAROpsToLs(NewOps, NewLengths, NewQL, NewTL);
		asserta(NewQL == OldQL);
		asserta(NewTL == OldTL);

		Ops = NewOps;
		Lengths = NewLengths;
		}

// 100M 6D M1
	if (Ops[N-1] == 'M' && Lengths[N-1] <= 2 && Lengths[N-2] > 4 && Ops[N-3] == 'M')
		{
		unsigned OldQL;
		unsigned OldTL;
		CIGAROpsToLs(Ops, Lengths, OldQL, OldTL);

		vector<char> NewOps;
		vector<unsigned> NewLengths;
		for (unsigned i = 0; i < N-1; ++i)
			{
			NewOps.push_back(Ops[i]);
			NewLengths.push_back(Lengths[i]);
			}
		NewLengths[N-3] += Lengths[N-1];

		unsigned NewQL;
		unsigned NewTL;
		CIGAROpsToLs(NewOps, NewLengths, NewQL, NewTL);
		asserta(NewQL == OldQL);
		asserta(NewTL == OldTL);

		Ops = NewOps;
		Lengths = NewLengths;
		}
	}

void OpsToCIGAR(const vector<char> &Ops, const vector<unsigned> &Lengths,
  string &CIGAR)
	{
	CIGAR.clear();
	const unsigned N = SIZE(Ops);
	asserta(SIZE(Lengths) == N);
	for (unsigned i = 0; i < N; ++i)
		Psa(CIGAR, "%u%c", Lengths[i], Ops[i]);
	}
