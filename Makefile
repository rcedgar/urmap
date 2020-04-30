CPP = ccache g++
CPPOPTS = -fopenmp -msse -mfpmath=sse -O3 -DNDEBUG -c

CC = ccache gcc
CCOPTS = -fopenmp -msse -mfpmath=sse -O3 -DNDEBUG -c

LNK = g++
LNKOPTS = -O3 -fopenmp -pthread -lpthread -static

HDRS = \
  aligner.h \
  alignresult.h \
  alnheuristics.h \
  alnparams.h \
  alpha.h \
  bitvec.h \
  chainer.h \
  chainer1.h \
  cigar.h \
  cmd.h \
  cmds.h \
  countsort.h \
  crc32.h \
  deflate.h \
  diagbox.h \
  estats.h \
  evalue.h \
  fastaseqsource.h \
  fastq.h \
  fastqseqsource.h \
  fileseqsource.h \
  filetype.h \
  getticks.h \
  globalaligner.h \
  gobuff.h \
  gzguts.h \
  hsp.h \
  hspfinder.h \
  inffast.h \
  inffixed.h \
  inflate.h \
  inftrees.h \
  label.h \
  linereader.h \
  localaligner.h \
  localaligner2.h \
  lockobj.h \
  lockobjs.h \
  mx.h \
  myopts.h \
  myutils.h \
  obj.h \
  objmgr.h \
  objtype.h \
  objtypes.h \
  omplock.h \
  outfiles.h \
  pathinfo.h \
  primes.h \
  readsimbench.h \
  samrec.h \
  samrec2.h \
  samseqsource.h \
  seqdb.h \
  seqdbseqsource.h \
  seqhash.h \
  seqinfo.h \
  seqsource.h \
  sort.h \
  sparsemx.h \
  state1.h \
  state2.h \
  svnversion.h \
  tenx.h \
  tracebit.h \
  trees.h \
  ufihit.h \
  ufihsp.h \
  ufindex.h \
  xdpmem.h \
  xtype.h \
  zconf.h \
  zlib.h \
  zutil.h \

OBJS = \
  o/alignresult.o \
  o/alnheuristics.o \
  o/alnparams.o \
  o/alpha.o \
  o/alpha2.o \
  o/arscorer.o \
  o/bitvec.o \
  o/chainer.o \
  o/chainer1.o \
  o/cigar.o \
  o/cmd.o \
  o/cmdline.o \
  o/diagbox.o \
  o/estats.o \
  o/evalue.o \
  o/fastaseqsource.o \
  o/fastq.o \
  o/fastqseqsource.o \
  o/filetype.o \
  o/fileseqsource.o \
  o/getcmd.o \
  o/gethsps.o \
  o/globalaligner.o \
  o/globalalignmem.o \
  o/gzipfileio.o \
  o/hspfinder.o \
  o/label.o \
  o/linereader.o \
  o/localaligner.o \
  o/localaligner2.o \
  o/localmulti.o \
  o/lockobj.o \
  o/logaln.o \
  o/makebitvec.o \
  o/mx.o \
  o/myutils.o \
  o/objmgr.o \
  o/onemasks.o \
  o/output2.o \
  o/outfiles.o \
  o/outputtab2.o \
  o/progress.o \
  o/rcemalloc.o \
  o/sam2aln.o \
  o/samrec2.o \
  o/samseqsource.o \
  o/searchbitvec.o \
  o/searchbitvec2.o \
  o/seqhash.o \
  o/setsam.o \
  o/pathinfo.o \
  o/pcb.o \
  o/prime.o \
  o/samrec.o \
  o/seqdb.o \
  o/seqdbfromfasta.o \
  o/seqdbio.o \
  o/seqdbseqsource.o \
  o/seqinfo.o \
  o/seqsource.o \
  o/setnucmx.o \
  o/tenx.o \
  o/test.o \
  o/timing.o \
  o/tracebackbitmem.o \
  o/map.o \
  o/alignhsp.o \
  o/alignhsp5.o \
  o/extendpen.o \
  o/extendscan.o \
  o/getseed.o \
  o/ufigetboth1s.o \
  o/map2.o \
  o/ufindex.o \
  o/ufindexio.o \
  o/output1.o \
  o/scan.o \
  o/scanslots.o \
  o/search1.o \
  o/search1m1.o \
  o/search1m6.o \
  o/state1.o \
  o/state2.o \
  o/search2.o \
  o/search2m4.o \
  o/search2m5.o \
  o/searchpending5.o \
  o/search1pepend.o \
  o/ufistats.o \
  o/viterbi.o \
  o/ungappedblast.o \
  o/urmap_main.o \
  o/viterbifastbandmem.o \
  o/viterbifastmem.o \
  o/xdropalignmem.o \
  o/xdropbwdmem.o \
  o/xdropbwdsplit.o \
  o/xdropfwdmem.o \
  o/xdropfwdsplit.o \
  o/adler32.o \
  o/crc32.o \
  o/deflate.o \
  o/gzlib.o \
  o/gzread.o \
  o/infback.o \
  o/inffast.o \
  o/inflate.o \
  o/inftrees.o \
  o/trees.o \
  o/zutil.o \

urmap : o/ $(OBJS)
	$(LNK) $(LNKOPTS) $(OBJS) -o o/urmap
	strip -d o/urmap

o/ :
	mkdir -p o/

o/adler32.o : adler32.c $(HDRS)
	$(CC) $(CCOPTS) -o o/adler32.o adler32.c

o/crc32.o : crc32.c $(HDRS)
	$(CC) $(CCOPTS) -o o/crc32.o crc32.c

o/deflate.o : deflate.c $(HDRS)
	$(CC) $(CCOPTS) -o o/deflate.o deflate.c

o/gzlib.o : gzlib.c $(HDRS)
	$(CC) $(CCOPTS) -o o/gzlib.o gzlib.c

o/gzread.o : gzread.c $(HDRS)
	$(CC) $(CCOPTS) -o o/gzread.o gzread.c

o/infback.o : infback.c $(HDRS)
	$(CC) $(CCOPTS) -o o/infback.o infback.c

o/inffast.o : inffast.c $(HDRS)
	$(CC) $(CCOPTS) -o o/inffast.o inffast.c

o/inflate.o : inflate.c $(HDRS)
	$(CC) $(CCOPTS) -o o/inflate.o inflate.c

o/inftrees.o : inftrees.c $(HDRS)
	$(CC) $(CCOPTS) -o o/inftrees.o inftrees.c

o/trees.o : trees.c $(HDRS)
	$(CC) $(CCOPTS) -o o/trees.o trees.c

o/zutil.o : zutil.c $(HDRS)
	$(CC) $(CCOPTS) -o o/zutil.o zutil.c

o/alignresult.o : alignresult.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/alignresult.o alignresult.cpp

o/alnheuristics.o : alnheuristics.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/alnheuristics.o alnheuristics.cpp

o/alnparams.o : alnparams.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/alnparams.o alnparams.cpp

o/alpha.o : alpha.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/alpha.o alpha.cpp

o/alpha2.o : alpha2.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/alpha2.o alpha2.cpp

o/arscorer.o : arscorer.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/arscorer.o arscorer.cpp

o/bitvec.o : bitvec.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/bitvec.o bitvec.cpp

o/chainer.o : chainer.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/chainer.o chainer.cpp

o/chainer1.o : chainer1.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/chainer1.o chainer1.cpp

o/cigar.o : cigar.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/cigar.o cigar.cpp

o/cmd.o : cmd.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/cmd.o cmd.cpp

o/cmdline.o : cmdline.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/cmdline.o cmdline.cpp

o/diagbox.o : diagbox.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/diagbox.o diagbox.cpp

o/estats.o : estats.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/estats.o estats.cpp

o/evalue.o : evalue.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/evalue.o evalue.cpp

o/fastaseqsource.o : fastaseqsource.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/fastaseqsource.o fastaseqsource.cpp

o/fastq.o : fastq.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/fastq.o fastq.cpp

o/fastqseqsource.o : fastqseqsource.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/fastqseqsource.o fastqseqsource.cpp

o/filetype.o : filetype.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/filetype.o filetype.cpp

o/fileseqsource.o : fileseqsource.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/fileseqsource.o fileseqsource.cpp

o/getcmd.o : getcmd.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/getcmd.o getcmd.cpp

o/gethsps.o : gethsps.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/gethsps.o gethsps.cpp

o/globalaligner.o : globalaligner.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/globalaligner.o globalaligner.cpp

o/globalalignmem.o : globalalignmem.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/globalalignmem.o globalalignmem.cpp

o/gzipfileio.o : gzipfileio.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/gzipfileio.o gzipfileio.cpp

o/hspfinder.o : hspfinder.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/hspfinder.o hspfinder.cpp

o/label.o : label.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/label.o label.cpp

o/linereader.o : linereader.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/linereader.o linereader.cpp

o/localaligner.o : localaligner.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/localaligner.o localaligner.cpp

o/localaligner2.o : localaligner2.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/localaligner2.o localaligner2.cpp

o/localmulti.o : localmulti.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/localmulti.o localmulti.cpp

o/lockobj.o : lockobj.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/lockobj.o lockobj.cpp

o/logaln.o : logaln.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/logaln.o logaln.cpp

o/makebitvec.o : makebitvec.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/makebitvec.o makebitvec.cpp

o/mx.o : mx.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/mx.o mx.cpp

o/myutils.o : myutils.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/myutils.o myutils.cpp

o/objmgr.o : objmgr.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/objmgr.o objmgr.cpp

o/onemasks.o : onemasks.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/onemasks.o onemasks.cpp

o/output2.o : output2.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/output2.o output2.cpp

o/outfiles.o : outfiles.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/outfiles.o outfiles.cpp

o/outputtab2.o : outputtab2.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/outputtab2.o outputtab2.cpp

o/progress.o : progress.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/progress.o progress.cpp

o/rcemalloc.o : rcemalloc.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/rcemalloc.o rcemalloc.cpp

o/sam2aln.o : sam2aln.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/sam2aln.o sam2aln.cpp

o/samrec2.o : samrec2.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/samrec2.o samrec2.cpp

o/samseqsource.o : samseqsource.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/samseqsource.o samseqsource.cpp

o/searchbitvec.o : searchbitvec.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/searchbitvec.o searchbitvec.cpp

o/searchbitvec2.o : searchbitvec2.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/searchbitvec2.o searchbitvec2.cpp

o/seqhash.o : seqhash.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/seqhash.o seqhash.cpp

o/setsam.o : setsam.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/setsam.o setsam.cpp

o/pathinfo.o : pathinfo.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/pathinfo.o pathinfo.cpp

o/pcb.o : pcb.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/pcb.o pcb.cpp

o/prime.o : prime.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/prime.o prime.cpp

o/samrec.o : samrec.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/samrec.o samrec.cpp

o/seqdb.o : seqdb.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/seqdb.o seqdb.cpp

o/seqdbfromfasta.o : seqdbfromfasta.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/seqdbfromfasta.o seqdbfromfasta.cpp

o/seqdbio.o : seqdbio.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/seqdbio.o seqdbio.cpp

o/seqdbseqsource.o : seqdbseqsource.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/seqdbseqsource.o seqdbseqsource.cpp

o/seqinfo.o : seqinfo.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/seqinfo.o seqinfo.cpp

o/seqsource.o : seqsource.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/seqsource.o seqsource.cpp

o/setnucmx.o : setnucmx.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/setnucmx.o setnucmx.cpp

o/tenx.o : tenx.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/tenx.o tenx.cpp

o/test.o : test.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/test.o test.cpp

o/timing.o : timing.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/timing.o timing.cpp

o/tracebackbitmem.o : tracebackbitmem.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/tracebackbitmem.o tracebackbitmem.cpp

o/map.o : map.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/map.o map.cpp

o/alignhsp.o : alignhsp.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/alignhsp.o alignhsp.cpp

o/alignhsp5.o : alignhsp5.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/alignhsp5.o alignhsp5.cpp

o/extendpen.o : extendpen.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/extendpen.o extendpen.cpp

o/extendscan.o : extendscan.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/extendscan.o extendscan.cpp

o/getseed.o : getseed.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/getseed.o getseed.cpp

o/ufigetboth1s.o : ufigetboth1s.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/ufigetboth1s.o ufigetboth1s.cpp

o/map2.o : map2.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/map2.o map2.cpp

o/ufindex.o : ufindex.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/ufindex.o ufindex.cpp

o/ufindexio.o : ufindexio.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/ufindexio.o ufindexio.cpp

o/output1.o : output1.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/output1.o output1.cpp

o/scan.o : scan.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/scan.o scan.cpp

o/scanslots.o : scanslots.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/scanslots.o scanslots.cpp

o/search1.o : search1.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/search1.o search1.cpp

o/search1m1.o : search1m1.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/search1m1.o search1m1.cpp

o/search1m6.o : search1m6.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/search1m6.o search1m6.cpp

o/state1.o : state1.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/state1.o state1.cpp

o/state2.o : state2.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/state2.o state2.cpp

o/search2.o : search2.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/search2.o search2.cpp

o/search2m4.o : search2m4.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/search2m4.o search2m4.cpp

o/search2m5.o : search2m5.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/search2m5.o search2m5.cpp

o/searchpending5.o : searchpending5.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/searchpending5.o searchpending5.cpp

o/search1pepend.o : search1pepend.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/search1pepend.o search1pepend.cpp

o/ufistats.o : ufistats.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/ufistats.o ufistats.cpp

o/viterbi.o : viterbi.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/viterbi.o viterbi.cpp

o/ungappedblast.o : ungappedblast.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/ungappedblast.o ungappedblast.cpp

o/urmap_main.o : urmap_main.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/urmap_main.o urmap_main.cpp

o/viterbifastbandmem.o : viterbifastbandmem.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/viterbifastbandmem.o viterbifastbandmem.cpp

o/viterbifastmem.o : viterbifastmem.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/viterbifastmem.o viterbifastmem.cpp

o/xdropalignmem.o : xdropalignmem.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/xdropalignmem.o xdropalignmem.cpp

o/xdropbwdmem.o : xdropbwdmem.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/xdropbwdmem.o xdropbwdmem.cpp

o/xdropbwdsplit.o : xdropbwdsplit.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/xdropbwdsplit.o xdropbwdsplit.cpp

o/xdropfwdmem.o : xdropfwdmem.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/xdropfwdmem.o xdropfwdmem.cpp

o/xdropfwdsplit.o : xdropfwdsplit.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/xdropfwdsplit.o xdropfwdsplit.cpp
