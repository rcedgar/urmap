#include "myutils.h"
#include "label.h"

void StripAllAnnots(string &Label)
	{
	size_t n = Label.find(';');
	if (n == string::npos || n == 0)
		return;
	Label = Label.substr(0, n);
	}

void GetAllAnnots(const string &Label, string &Annots)
	{
	Annots.clear();
	size_t n = Label.find(';');
	if (n == string::npos || n == 0)
		return;
	unsigned L = SIZE(Label);
	for (unsigned i = unsigned(n) + 1; i < L; ++i)
		{
		char c = Label[i];
		Annots.push_back(c);
		}
	}

const char *GetStrField(const string &Label, const string &NameEq,
  string &Value)
	{
	Value.clear();
	vector<string> Fields;
	Split(Label, Fields, ';');
	const unsigned N = SIZE(Fields);
	for (unsigned i = 0; i < N; ++i)
		{
		const string &Field = Fields[i];
		if (StartsWith(Field, NameEq))
			{
			Value = Field.substr(NameEq.size(), string::npos);
			break;
			}
		}
	return Value.c_str();
	}

void StripAnnot(string &Label, const string &NameEq)
	{
	if (Label.find(NameEq) == string::npos)
		return;

	string NewLabel;
	vector<string> Fields;
	Split(Label, Fields, ';');
	const unsigned N = SIZE(Fields);
	for (unsigned i = 0; i < N; ++i)
		{
		const string &Field = Fields[i];
		if (Field.find(NameEq) == 0)
			continue;
		NewLabel += Field + ";";
		}
	if (NewLabel.find('=') == string::npos)
		{
		Label.clear();
		unsigned n = SIZE(NewLabel);
		for (unsigned i = 0; i + 1 < n; ++i)
			Label.push_back(NewLabel[i]);
		}
	else
		Label = NewLabel;
	}

void AppendIntField(string &Label, const string &NameEq, unsigned Value)
	{
	Psasc(Label, "%s%u", NameEq.c_str(), Value);
	}

void AppendStrField(string &Label, const string &NameEqValue)
	{
	Psasc(Label, "%s", NameEqValue.c_str());
	}

void AppendStrField(string &Label, const string &NameEq, const string &Value)
	{
	Psasc(Label, "%s%s", NameEq.c_str(), Value.c_str());
	}

void AppendSize(string &Label, unsigned Size)
	{
	AppendIntField(Label, "size=", Size);
	}

void AppendTaxStr(string &Label, const string &s)
	{
	AppendStrField(Label, "tax=", s);
	}

void StripSize(string &Label)
	{
	StripAnnot(Label, "size=");
	}

void StripTax(string &Label)
	{
	StripAnnot(Label, "tax=");
	StripAnnot(Label, "taxstr=");
	StripAnnot(Label, "sciname=");
	StripAnnot(Label, "taxlev=");
	size_t n = Label.find("\tRoot;");
	if (n != string::npos)
		Label = Label.substr(0, n);
	}

const char *GetTaxStrFromLabel(const string &Label, string &s)
	{
	GetStrField(Label, "tax=", s);
	return s.c_str();
	}

void ReplaceSizeInLabel(string &Label, unsigned NewSize)
	{
	StripSize(Label);
	AppendSize(Label, NewSize);
	}

unsigned GetIntFieldFromLabel(const string &Label, const string &NameEq, unsigned Default)
	{
	const unsigned n = SIZE(NameEq);
	vector<string> Fields;
	Split(Label, Fields, ';');
	const unsigned N = SIZE(Fields);
	for (unsigned i = 0; i < N; ++i)
		{
		const string &Field = Fields[i];
		if (Field.substr(0, n) == NameEq)
			{
			string a = Field.substr(n, string::npos);
			const char *as = a.c_str();
			if (!IsUintStr(as) && Default == UINT_MAX)
				Die("%s not integer >%s", NameEq.c_str(), Label.c_str());
			if (!IsUintStr(as))
				Die("%s invalid value >%s", NameEq.c_str(), Label.c_str());
			unsigned Value = StrToUint(as);
			return Value;
			}
		}
	if (Default == UINT_MAX)
		Die("%s not found in label >%s", NameEq.c_str(), Label.c_str());
	return Default;
	}

unsigned GetSizeFromLabel(const string &Label, unsigned Default)
	{
	StartTimer(GetSizeFromLabel);
//	unsigned Size = GetIntFieldFromLabel(Label, "size=", Default);
	unsigned Size = Default;
	const char *p = strstr(Label.c_str(), ";size=");
	if (p != 0)
		Size = (unsigned) atoi(p+6);
	else if (Default == UINT_MAX)
		Die("Missing size= in >%s", Label.c_str());
	EndTimer(GetSizeFromLabel);
	return Size;
	}

unsigned GetTaxIdFromLabel(const string &Label, unsigned Default)
	{
	return GetIntFieldFromLabel(Label, "tax=", Default);
	}

void GetAccFromLabel(const string &Label, string &Acc)
	{
	Acc.clear();
	const unsigned N = SIZE(Label);
	for (unsigned i = 0; i < N; ++i)
		{
		char c = Label[i];
		if (c == ' ' || c == '|' || c == ';')
			{
			if (Acc != "gi")
				return;
			}
		Acc += c;
		}
	}

byte GetSplitRankFromLabel(const string &Label)
	{
	string s;
	GetStrField(Label, "split=", s);
	if (SIZE(s) > 1)
		Die("Invalid split in label >%s", Label.c_str());
	return s.c_str()[0];
	}

void GetOTUNameFromLabel(const string &Label, string &OTUName)
	{
	GetStrField(Label, "otu=", OTUName);
	if (!OTUName.empty())
		return;
	
	GetAccFromLabel(Label, OTUName);
	if (OTUName.empty())
		Die("Empty OTU name in label >%s", Label.c_str());
	}

void GetSampleNameFromLabel(const string &Label, string &SampleName)
	{
	SampleName.clear();

	GetStrField(Label, "sample=", SampleName);
	if (!SampleName.empty())
		return;

	GetStrField(Label, "barcodelabel=", SampleName);
	if (!SampleName.empty())
		return;

	if (optset_sample_delim)
		{
		const string &d = opt(sample_delim);
		size_t n = Label.find(d);
		if (n == string::npos)
			Die("delim '%s' not found in >%s", d.c_str(), Label.c_str());
		SampleName = Label.substr(0, n);
		return;
		}

	unsigned L = SIZE(Label);
	for (unsigned i = 0; i < L; ++i)
		{
		char c = Label[i];
		if (!isalpha(c) && !isdigit(c) && !(c == '_'))
			return;
		SampleName.push_back(c);
		}
	}

// urtils/varsim.cpp:142
//              0  1  2  3  4  5  6  7             0     1        2            3
// Ps(Label, "R%u.%u.%u.%u.%s.%c.%u.%u", m_PairIndex, m_BC, m_MolId, PairIndex+1,
//   m_RefLabel.c_str(), m_AorB, R1Pos+1, R2Pos+1);
//            4               5        6        7
void ParseVarsimLabel(const string &ReadLabel, unsigned &ReadPairIndex,
  uint32 &BC, unsigned &MolId, unsigned &MolPairIndex,
  string &ChrName, char &AorB, unsigned &R1Pos, unsigned &R2Pos)
	{
	string Label = ReadLabel;
	TruncWhiteSpace(Label);
	vector<string> Fields;
	Split(Label, Fields, '.');
	asserta(SIZE(Fields) == 8);
	asserta(Label[0] == 'R');
	ReadPairIndex = StrToUint(Fields[0].c_str() + 1);
	BC = StrToUint(Fields[1]);
	MolId = StrToUint(Fields[2]);
	MolPairIndex = StrToUint(Fields[3]);
	ChrName = Fields[4];
	string sAorB = Fields[5];
	AorB = sAorB[0];
	if (AorB != 'A' && AorB != 'B')
		Die("Invalid label, AorB %s", ReadLabel.c_str());
	R1Pos = StrToUint(Fields[6]);
	R2Pos = StrToUint(Fields[7]);
	}

uint32 GetBCFromLabel(const string &ReadLabel)
	{
	unsigned ReadPairIndex, MolId, MolPairIndex, R1Pos, R2Pos;
	char AorB;
	unsigned n = SIZE(ReadLabel);
	//asserta(ReadLabel[n-2] == ' ');
	//char OneTwo = ReadLabel[n-1];
	string TargetLabel;
	uint32 BC;
	ParseVarsimLabel(ReadLabel, ReadPairIndex, BC, MolId, MolPairIndex,
	  TargetLabel, AorB, R1Pos, R2Pos);
	return BC;
	}

void GetCoordsFromLabel_varsim(const string &ReadLabel,
  string &TargetLabel, unsigned &Pos)
	{
	unsigned ReadPairIndex, MolId, MolPairIndex, R1Pos, R2Pos;
	char AorB;
	unsigned n = SIZE(ReadLabel);
	asserta(ReadLabel[n-2] == ' ');
	char OneTwo = ReadLabel[n-1];
	uint32 BC;
	ParseVarsimLabel(ReadLabel, ReadPairIndex, BC, MolId, MolPairIndex,
	  TargetLabel, AorB, R1Pos, R2Pos);
	if (OneTwo == '1')
		Pos = R1Pos;
	else if (OneTwo == '2')
		Pos = R2Pos;
	else
		Die("OneTwo %c", OneTwo);
	}

// Sim.5;label=chr1;strand=-;pos=3411;diffs=S28TC;
// Pos is 1-based lo coord w.r.t. plus strand of reference
void GetCoordsFromLabel(const string &Label,
  string &TargetLabel, unsigned &Pos, bool &Strand)
	{
	TargetLabel.clear();
	Pos = UINT_MAX;
	vector<string> Fields;
	Split(Label, Fields, ';');
	const unsigned n = SIZE(Fields);
	vector<string> Fields2;
	for (unsigned i = 1; i < n; ++i)
		{
		const string &Fieldi = Fields[i];
		Split(Fieldi, Fields2, '=');
		if (SIZE(Fields2) != 2)
			continue;
		const string &Name = Fields2[0];
		const string &Value = Fields2[1];
		if (Name == "label")
			TargetLabel = Value;
		else if (Name == "pos")
			Pos = StrToUint(Value);
		else if (Name == "strand")
			{
			if (Value == "+")
				Strand = true;
			else if (Value == "-")
				Strand = false;
			else
				asserta(false);
			}
		}
	if (Pos == UINT_MAX)
		Die("Pos=UINT_MAX, label >%s", Label.c_str());
	}
