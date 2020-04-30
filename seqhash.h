#ifndef seqhash_h
#define seqhash_h

class SeqDB;

uint32 SeqHash32(const byte *Seq, unsigned L);
uint32 SeqHashRC32(const byte *Seq, unsigned L);
uint32 SeqHash32_EitherStrand(const byte *Seq, unsigned L);
bool SeqEq(const byte *Seq1, unsigned L1, const byte *Seq2, unsigned L2);
bool SeqEqRC(const byte *Seq1, unsigned L1, const byte *Seq2, unsigned L2);

#endif // seqhash_h
