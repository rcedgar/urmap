#ifndef readsimbench_h
#define readsimbench_h

struct RSBInfo
	{
	unsigned V;
	char ab;
	unsigned Index;
	string TargetLabel;
	bool Plus;
	unsigned TargetLo;
	unsigned VPos;
	bool IsFwd;
	bool HasV;
	string Aln;
	string SubErrs;

	RSBInfo()
		{
		Clear();
		}

	void Clear()
		{
		V = UINT_MAX;
		ab = '?';
		Index = UINT_MAX;
		TargetLabel.clear();
		Plus = false;
		TargetLo = UINT_MAX;
		VPos = UINT_MAX;
		IsFwd = false;
		HasV = false;
		Aln.clear();
		SubErrs.clear();
		}

	void LogMe() const
		{
		Log("RSB: V=%u", V);
		Log(", ab=%c", ab);
		Log(", Index=%u", Index);
		Log(", target=%s(%u)%c", TargetLabel.c_str(), TargetLo+1, pom(Plus));
		Log(", %s", IsFwd ? "Fwd" : "Rev");
		Log(", hasv=%c", yon(HasV));
		Log(", aln=%s", Aln.c_str());
		Log(", suberrs=%s", SubErrs.c_str());
		Log("\n");
		}
	};

void GetCoordsFromLabel(const string &Label,
  string &TargetLabel, unsigned &Pos, bool &Strand);

void GetVarIxFromLabel(const string &Label, string &VarIx);

void GetCoordsFromLabel_varsim(const string &ReadLabel,
  string &TargetLabel, unsigned &Pos);

uint32 GetBCFromLabel(const string &ReadLabel);

void ParseVarsimLabel(const string &ReadLabel, unsigned &ReadPairIndex,
  uint32 &BC, unsigned &MolId, unsigned &MolPairIndex,
  string &ChrName, char &AorB, unsigned &R1Pos, unsigned &R2Pos);

#endif // readsimbench_h
