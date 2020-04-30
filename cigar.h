#ifndef cigar_h
#define cigar_h

void PathToCIGAR(const char *Path, string &CIGAR);
unsigned CIGARToQL(const string &CIGAR);
void CIGARToLs(const string &CIGAR, unsigned &QL, unsigned &TL);
void CIGAROpsToLs(const vector<char> &Ops, const vector<unsigned> &Lengths,
  unsigned &QL, unsigned &TL);
void CIGARGetOps(const string &CIGAR, vector<char> &Ops,
  vector<unsigned> &Lengths);
void OpsToCIGAR(const vector<char> &Ops, const vector<unsigned> &Lengths,
  string &CIGAR);
void CIGAROpsFixDanglingMs(vector<char> &Ops, vector<unsigned> &Lengths);

#endif // cigar_h
