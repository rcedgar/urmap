#ifndef tenx_h
#define tenx_h

#include <set>

typedef uint32 count_t;
const uint32 TX_MAXCOUNT = uint32(count_t(~0));
const uint64 TX_NBC = (uint64(1) << 32);

/***
https://support.10xgenomics.com/genome-exome/software/pipelines/latest/output/bam

	Tag	 Desc
	 BX  error-corrected barcode
	 RX  raw barcode
	 QX  barcode qual
	 TR  trimmed 7mer
	 TQ  trimmed qual
***/

uint32 Tenx_GetIntegerBarcode(const byte *Seq);
const char *Tenx_GetBarcodeFromSeq(const byte *Seq, string &BC);
uint32 Tenx_GetIntegerBarcodeFromLabel(const char *Label);
void Tenx_ReadWhitelist(const string &FileName, set<uint32> &WhiteSet);
const char *Tenx_IntegerBarcodeToStr(uint32 BC, string &s);
void Tenx_ClusterCoords(uint32 BC, const vector <uint32> &Coords);

#endif // tenx_h
