/*  File: synUtil.c
 *  Author: Chenxi Zhou (cz370@cam.ac.uk)
 *  Copyright (C) Chenxi Zhou, Cambridge University, 2023
 *-------------------------------------------------------------------
 * Description: syncmer utilities
 * Exported functions:
 * HISTORY:
 * Last edited: Oct 15 22:59 2024 (cz370)
 * Created: Tue Oct 12 11:57:13 2024 (cz370)
 *-------------------------------------------------------------------
 */

#define SYNC_VERSION "0.0"

#include <unistd.h>

#include "ONElib.h"
#include "seqhash.h"
#include "syng.h"
#include "cn_plot.R.h"

static void setParams (OneFile *of)
{ // set kmer size params
  params.k = oneInt(of, 0) ;
  params.w = oneInt(of, 1) ;
  params.seed = oneInt(of,2) ;
}

static void checkParams (OneFile *of)
{ // check consistency
  if (oneInt(of,0) != params.k ||
      oneInt(of,1) != params.w ||
      oneInt(of,2) != params.seed)
    die ("hash parameters mismatch: (k,w,s) file (%d,%d,%d) != code (%d,%d,%d)",
        oneInt(of,0), oneInt(of,1), oneInt(of,2), params.k, params.w, params.seed) ;
}

void SystemX(char *command)
{ // system call
  if (system(command) != 0)
    die ("error running system command: %s",command) ;
}

int main_fish  (int argc, char *argv[]) ;
int main_hist  (int argc, char *argv[]) ;
int main_hist2 (int argc, char *argv[]) ;

static char usage[] =
  "Usage: synUtil <options>\n"
  "  fish     :  extract sequences by syncmer fishing\n"
  "  hist     :  make syncmer count histogram from a khash file\n"
  "  hist2    :  make two-way syncmer count histogram from two khash files\n"
  ;

static char hist_usage[] =
  "Usage: synUtil hist [options] seq.1khash\n"
  "  -c <max syncmer copy>    : 0\n"
  "  -o <output file>         : defualt to stdout\n"
  ;

static char hist2_usage[] =
  "Usage: synUtil hist2 [options] seq1.1khash seq2.1khash\n"
  "  -c1 <max syncmer copy>   : max syncmer copy in seq1 [6]\n"
  "  -c2 <max syncmer copy>   : max syncmer copy in seq2 [300]\n"
  "  -p                       : make a line plot (need -o)\n"
  "  -v                       : verbose mode\n"
  "  -o  <output file>        : defualt to stdout\n"
  ;

static char fish_usage[] =
  "Usage: synUtil fish [options] ref.1khash seq.1khash seq.1syncseq\n"
  "  -c1 <min syncmer copy>   : min syncmer copy in sequence [3]\n"
  "  -c2 <max syncmer copy>   : max syncmer copy in sequence\n"
  "  -a  <min average>        : min average syncmer copy [0]\n"
  "  -f  <max low-copy frac>  : max low-copy syncmer fraction [1.0]\n"
  "  -n  <min syncmer number> : min syncmer hits in sequence [1]\n"
  "  -o  <output file>        : defualt to stdout\n"
  ;

int main(int argc, char *argv[])
{
  argc-- ; ++argv;
  if (!argc) { printf ("%s",usage) ; exit (0) ; }
  if (strcmp(*argv, "fish") == 0) return main_fish(argc, argv) ;
  if (strcmp(*argv, "hist") == 0) return main_hist(argc, argv) ;
  if (strcmp(*argv, "hist2") == 0) return main_hist2(argc, argv) ;
  else {
		die ("unrecognized command '%s'\n%s", *argv, usage) ;
		return 1 ;
	}
}

typedef struct {
  I64 index ;
  I64 count ; 
} KmerCount ;

static int KSORT(const void *a, const void *b)
{ I64 x, y;
  x = ((KmerCount *) a)->index;
  y = ((KmerCount *) b)->index;
  return (x > y) - (x < y) ;
}

int main_hist (int argc, char *argv[])
{
  char *outFileName, *syncFileName ;
  int maxC ;
  FILE *out ;
  Array kmerCounts ;
  Hash kmerHash ;
  
  maxC = 0 ;
  kmerCounts = NULL ;
  kmerHash = NULL;
  outFileName = NULL ;
  syncFileName = NULL ;
  
  timeUpdate (0) ;

  storeCommandLine (argc, argv) ;
  argc-- ; ++argv ;
  if (!argc) { printf ("%s",hist_usage) ; exit (0) ; }
  while (argc > 0 && **argv == '-')
    if (!strcmp (*argv, "-c") && argc > 1) { maxC = atoi(argv[1]) ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (*argv, "-o") && argc > 1) { outFileName = argv[1] ; argc -= 2 ; argv += 2 ; }
    else die ("unknown parameter %s\n%s", *argv, hist_usage) ;
  if (!maxC)
    { // use a table for the syncmers of counts < 65536
      // and a hash table for >= 65536
      maxC = 1<<16 ;
      kmerCounts = arrayCreate (1<<16, KmerCount) ;
      kmerHash = hashCreate (1<<16) ;
    }
  if (argc < 1)
    die ("need an input syncmer hash table\n%s", hist_usage) ;
  syncFileName = argv[0] ;

  U64 totSync = 0 ;
  OneSchema  *schema = oneSchemaCreateFromText (syngSchemaText) ;
  
  OneFile *of = oneFileOpenRead (syncFileName, schema, "khash", 1) ;
  if (!of) die ("failed to open syncmer file %s to read", syncFileName) ;
  
  // read the syncmer hash parameters
  while (oneReadLine (of) && of->lineType != 'h') ;
  if (of->lineType == 'h') setParams (of) ;
  else die ("sync file %s has no 'h' parameters record", syncFileName) ;

  // get 't' line
  while (oneReadLine (of) && of->lineType != 't') ;
  if (of->lineType != 't') die ("failed to find 't' line in khash file") ;
  totSync = oneInt(of,0) ;

  // read in the node counts 
  while (oneReadLine (of) && of->lineType != 'C') ;
  if (of->lineType != 'C' || oneLen(of) != totSync)
	  die ("failed to find expected C line in khash file") ;
  I64 i, c, *hist = new0 (maxC+1, I64), *counts = oneIntList(of) ;
  int index;
  for (i = 0 ; i < totSync ; ++i) {
    c = counts[i] ;
    if (c <= maxC)
      hist[c-1] += 1 ;
    else if (kmerHash)
      { if (hashAdd (kmerHash, HASH_INT(c), &index))
          array (kmerCounts, index, KmerCount) = (KmerCount) {c, 1} ;
        else
          arrayp (kmerCounts, index, KmerCount)->count += 1 ;
      }
    else
      hist[maxC] += 1 ;
  }

  if (outFileName)
    { out = fopen(outFileName, "w") ;
      if (!out) die("cannot open file '%s' to write\n", outFileName) ;
    }
  else out = stdout ;
  for (i = 0; i < maxC; i++)
    if (hist[i])
      fprintf(out, "%12lld %12lld\n", i+1, hist[i]) ;
  if (kmerHash)
    { c = hashCount (kmerHash) ;
      qsort (kmerCounts->base, c, sizeof(KmerCount), KSORT) ;
      for (i = 0; i < c; i++)
        fprintf(out, "%12lld %12lld\n", 
          arrayp (kmerCounts, i, KmerCount)->index,
          arrayp (kmerCounts, i, KmerCount)->count) ;
    }
  else if (hist[maxC])
      fprintf(out, "%12d %12lld\n", maxC+1, hist[maxC]) ;

  if (out != stdout)
    fclose (out) ;
  oneSchemaDestroy (schema) ;
  oneFileClose (of) ; // close the old file

  free (hist) ;
  if (kmerCounts)
    arrayDestroy (kmerCounts) ;
  if (kmerHash)
    hashDestroy (kmerHash) ;
  
  free(getCommandLine()) ;

  timeUpdate (stderr) ;
  fprintf (stderr, "total: ") ; timeTotal (stderr) ;
  return 0 ;
}

#define min(x,y) ((x)<(y)? (x) : (y))

int main_hist2 (int argc, char *argv[])
{
  char *outFileName, *aSyncFileName, *rSyncFileName ;
  int maxA, maxR ;
  FILE *outCNI ;
  bool plot, verb ;
  
  maxA = 6 ;
  maxR = 300 ;
  plot = false ;
  verb = false ;
  outFileName = NULL ;

  timeUpdate (0) ;
  
  storeCommandLine (argc, argv) ;
  argc-- ; ++argv ;
  if (!argc) { printf ("%s",hist2_usage) ; exit (0) ; }
  while (argc > 0 && **argv == '-')
    if (!strcmp (*argv, "-c1") && argc > 1) { maxA = atoi(argv[1]) ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (*argv, "-c2") && argc > 1) { maxR = atoi(argv[1]) ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (*argv, "-p")) { plot = true ; argc -= 1 ; argv += 1 ; }
    else if (!strcmp (*argv, "-v")) { verb = true ; argc -= 1 ; argv += 1 ; }
    else if (!strcmp (*argv, "-o") && argc > 1) { outFileName = argv[1] ; argc -= 2 ; argv += 2 ; }
    else die ("unknown parameter %s\n%s", *argv, hist2_usage) ;
  if (argc < 2)
    die ("need two input syncmer tables") ;
  aSyncFileName = argv[0] ;
  rSyncFileName = argv[1] ;
  if (!outFileName && plot) die ("need -o for plotting") ;

  OneSchema *schema ;
  OneFile *of ;
  int   i, j ;
  I64   iSync, nSync ;
  I64  *aSyncTable, **mSyncTable, *cref, *cseq ;
  KmerHash   *khref ;

  aSyncTable = new0((I64)(maxA+2)*(maxR+2), I64) ;
  mSyncTable = new0(maxA+2, I64 *) ;
  mSyncTable[0] = aSyncTable ;
  
  for (i = 1; i <= maxA+1; i++)
    mSyncTable[i] = mSyncTable[i-1] + maxR + 2 ;

  schema = oneSchemaCreateFromText (syngSchemaText) ;
  of = oneFileOpenRead (aSyncFileName, schema, "khash", 1) ;
  if (!of) die ("failed to open syncmer file %s to read", aSyncFileName) ;
  
  // read the syncmer hash parameters
  while (oneReadLine (of) && of->lineType != 'h') ;
  if (of->lineType == 'h') setParams (of) ;
  else die ("sync file %s has no 'h' parameters record", aSyncFileName) ;

  // read in the kmerhash
  khref = kmerHashReadOneFile (of) ;
  if (khref->len != params.w + params.k)
	  die ("syncmer len mismatch %d != %d + %d", khref->len, params.w, params.k) ;

  // read in the node counts 
  if (of->lineType != 'C' || oneLen(of) != kmerHashMax(khref))
    die ("failed to find expected C line in khash file %s", aSyncFileName) ;
  cref = new (khref->max, I64) ;
  memcpy (cref, oneIntList(of), sizeof(I64) * khref->max) ;

  if (verb) fprintf (stderr, "read %lld nodes from %s\n", kmerHashMax(khref), aSyncFileName) ;
  oneFileClose (of) ; // close the khash file

  // read in the read khash
  of = oneFileOpenRead (rSyncFileName, schema, "khash", 1) ;
  if (!of) die ("failed to open syncmer file %s to read", rSyncFileName) ;
  
  // consistency check for kmer size
  while (oneReadLine (of) && of->lineType != 'h') ;
  if (of->lineType == 'h') checkParams (of) ;
  else die ("sync file %s has no 'h' parameters record", rSyncFileName) ;

  if (!oneGoto (of, 'S', 0)) die ("failed to read 'S' line in khash file %s", rSyncFileName) ;
  while (of->lineType != 't' && oneReadLine (of)) ;
  if (of->lineType != 't') die ("failed to find 't' line in khash file %s", rSyncFileName) ;
  if (khref->plen != (oneInt(of,1)+31) >> 5) die ("different plen in two khash files") ;
  nSync = oneInt(of, 0) ;
  
  // go to the 'C' line for node counts
  while (of->lineType != 'C' && oneReadLine (of)) ;
  if (of->lineType != 'C' || oneLen(of) != nSync)
    die ("failed to find expected C line in khash file %s", rSyncFileName) ;
  cseq = new (nSync, I64) ;
  memcpy (cseq, oneIntList(of), sizeof(I64) * nSync) ;

  // read the DNA, maybe in chunks - assume they are at least multiple of 4bp, i.e. full bytes
  U8   *dna, *dna0, *dna1 ;
  I64   acnt, rcnt, index, resid ;
  int   plen ;
  bool  exist ;
  plen = khref->plen<<3 ; // in bytes instead of U64
  dna  = new(min(nSync*plen, ((I64)1<<27)+plen), U8) ; // extra plen for residuals
  iSync = 0 ;
  resid = 0 ;
  if (!oneGoto (of, 'S', 0)) die ("failed to read 'S' line in khash file %s", rSyncFileName) ;
  while (of->lineType != 't' && oneReadLine (of)) ;
  while (oneReadLine (of) && of->lineType == 'S')
    { memcpy(dna + resid, oneDNA2bit(of), oneLen(of)>>2) ;
      dna1 = dna + resid + (oneLen(of)>>2) - plen ;
      for (dna0 = dna; dna0 <= dna1; dna0 += plen)
        { rcnt = cseq[iSync] ;
          exist = kmerHashFindPacked (khref, (U64 *) dna0, &index) ;
          acnt = exist? cref[index-1] : 0 ;
          if (rcnt > maxR) rcnt = maxR + 1 ;
          if (acnt > maxA) acnt = maxA + 1 ;
          mSyncTable[acnt][rcnt] += 1 ;
          iSync++;
        }
      resid = dna1 + plen - dna0 ;
      memmove(dna, dna0, resid) ; // copy residuals
    }
  if (iSync != nSync) die ("wrong number of syncmer in kmerhash %lld %lld", iSync, nSync) ;
  free (dna) ;

  if (verb) fprintf (stderr, "read %lld nodes from %s\n", nSync, rSyncFileName) ;
  oneFileClose (of) ; // close the khash file

  // count asm syncs
  iSync = kmerHashMax(khref) ;
  while (iSync--)
    { acnt = cref[iSync] ;
      if (acnt > maxA) acnt = maxA + 1 ;
      mSyncTable[acnt][0] += 1 ;
    }
  for (i = 1; i <= maxA + 1; i++) 
    { aSyncTable = mSyncTable[i] ;
      acnt = 0;
      for (j = 1; j <= maxR + 1; j++)
        acnt += aSyncTable[j] ;
      assert(aSyncTable[0] >= acnt) ;
      aSyncTable[0] -= acnt ;
    }

  if (verb)
    for (i = 0; i <= maxR + 1; i++) 
      { fprintf(stderr, "%-4d:", i) ;
        for (j = 0; j <= maxA + 1; j++)
            fprintf(stderr, " %9lld", mSyncTable[j][i]) ;
        fprintf(stderr, "\n") ;
      }

  outCNI = NULL ;
  if (outFileName)
    { char *f = new (strlen(outFileName) + 6, char) ;
      sprintf(f, "%s.cni", outFileName) ;
      outCNI = fopen(f, "w") ;
      free(f) ;
    }
  else
    outCNI = stdout;
  fprintf(outCNI, "Copies\tkmer_multiplicity\tCount\n") ;
  for (j = 1; j <= maxR + 1; j++)
    fprintf(outCNI, "read-only\t%-9d\t%lld\n", j, mSyncTable[0][j]) ;
  for (i = 1; i <= maxA; i++)
    for (j = 1; j <= maxR + 1; j++)
      fprintf(outCNI, "%-9d\t%-9d\t%lld\n", i, j, mSyncTable[i][j]) ;
  for (j = 1; j <= maxR + 1; j++)
    fprintf(outCNI, ">%-8d\t%-9d\t%lld\n", maxA, j, mSyncTable[maxA+1][j]) ;
  if (outCNI != stdout)
    fclose(outCNI) ;
  
  if (plot)
    { double XDIM, YDIM ;
      int XMAX ;
      I64 YMAX ;
      
      YMAX = 0;
      for (i = 1; i <= maxA + 1; i++)
        for (j = 1; j <= maxR + 1; j++)
          if (mSyncTable[i][j] > YMAX)
            YMAX = mSyncTable[i][j] ;
      YMAX *= 1.4 ;

      nSync = kmerHashMax(khref) * .99 ;
      for (XMAX = 1; XMAX <= maxR + 1; XMAX++)
        { for (j = 1; j <= maxA + 1; j++)
            nSync -= mSyncTable[j][XMAX] ;
          if (nSync <= 0) break ;
        }
      
      XDIM = 8. ;
      YDIM = 6. ;

      FILE *outSRC ;
      int desc ;
      char *temp = new (strlen(outFileName) + 12, char) ;
      sprintf(temp, "%sXXXXXX", outFileName) ;
      desc = mkstemp(temp) ;
      if (desc == -1) die ("error creating temporary file") ;
      outSRC = fdopen(desc, "w") ;
      if (!outSRC) die ("error opening temporary file") ;
      fwrite(cn_plot_script, strlen(cn_plot_script), 1, outSRC) ;
      fclose(outSRC) ;

      char *cmd = new (strlen(outFileName)*3+256, char) ;
      sprintf(cmd,"Rscript %s -f %s.cni -o %s -p -x %g -y %g -m %d -n %lld -t line 2>/dev/null",
        temp,outFileName,outFileName,XDIM,YDIM,XMAX,YMAX) ;
      if (verb) fprintf(stderr, "run system command: %s\n", cmd) ;
      SystemX(cmd) ;
      if (remove(temp)) warn ("error deleting temporary file") ;
      
      free (temp) ;
      free (cmd) ;
    }

  kmerHashDestroy (khref) ;

  free (mSyncTable[0]) ;
  free (mSyncTable) ;
  free (cref) ;
  free (cseq) ;

  oneSchemaDestroy (schema) ;

  free(getCommandLine()) ;

  timeUpdate (stderr) ;
  fprintf (stderr, "total: ") ; timeTotal (stderr) ;
  return 0 ;
}

#define packseq(kh,i) ((kh)->pack + (i)*(kh)->plen)

#define PRINT_SYNC_PROFILE

int main_fish (int argc, char *argv[])
{
  char *outFileName, *syncRefHashFileName, *syncSeqHashFileName, *syncSeqFileName ;
  int minC, maxC, minN, minA ;
  double maxF ;
  FILE *out ;
    
  minC = 3 ;
  maxC = INT_MAX ;
  minN = 1 ;
  minA = 0 ;
  maxF = 1.0 ;
  outFileName = NULL ;
  syncRefHashFileName = NULL ;
  syncSeqHashFileName = NULL ;
  syncSeqFileName = NULL ;

  timeUpdate (0) ;

  storeCommandLine (argc, argv) ;
  argc-- ; ++argv ;
  if (!argc) { printf ("%s",fish_usage) ; exit (0) ; }
  while (argc > 0 && **argv == '-')
    if (!strcmp (*argv, "-c1") && argc > 1) { minC = atoi(argv[1]) ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (*argv, "-c2") && argc > 1) { maxC = atoi(argv[1]) ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (*argv, "-a") && argc > 1) { minA = atoi(argv[1]) ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (*argv, "-f") && argc > 1) { maxF = atof(argv[1]) ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (*argv, "-n") && argc > 1) { minN = atoi(argv[1]) ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (*argv, "-o") && argc > 1) { outFileName = argv[1] ; argc -= 2 ; argv += 2 ; }
    else die ("unknown parameter %s\n%s", *argv, fish_usage) ;
  if (!minC) minC = 3 ;
  if (!maxC) maxC = INT_MAX ;
  if (!minN) minN = 1 ;
  if (outFileName)
    { out = fopen(outFileName, "w") ;
      if (!out) die("cannot open file '%s' to write\n", outFileName) ;
    } 
  else out = stdout ;
  if (minC > maxC)
    die("min syncmer count greater than max") ;
  if (argc < 3)
    die ("need at least two input files\n%s", fish_usage) ;
  syncRefHashFileName = argv[0] ;
  syncSeqHashFileName = argv[1] ;
  syncSeqFileName = argv[2] ;

  U64 nSeq = 0, totSync = 0 ;
  I64        *cref, *cseq ;
  OneSchema  *schema = oneSchemaCreateFromText (syngSchemaText) ;
  KmerHash   *khref, *khseq ;

  OneFile *of = oneFileOpenRead (syncRefHashFileName, schema, "khash", 1) ;
  if (!of) die ("failed to open syncmer file %s to read", syncRefHashFileName) ;
  
  // read the syncmer hash parameters
  while (oneReadLine (of) && of->lineType != 'h') ;
  if (of->lineType == 'h') setParams (of) ;
  else die ("sync file %s has no 'h' parameters record", syncRefHashFileName) ;

  // read in the kmerhash
  khref = kmerHashReadOneFile (of) ;
  if (khref->len != params.w + params.k)
	  die ("syncmer len mismatch %d != %d + %d", khref->len, params.w, params.k) ;

  // read in the node counts 
  if (of->lineType != 'C' || oneLen(of) != kmerHashMax(khref))
    die ("failed to find expected C line in khash file") ;
  cref = new (khref->max, I64) ;
  memcpy (cref, oneIntList(of), sizeof(I64) * khref->max) ;

  fprintf (stderr, "read %llu nodes from %s\n", kmerHashMax(khref), syncRefHashFileName) ;
  oneFileClose (of) ; // close the khash file

  of = oneFileOpenRead (syncSeqHashFileName, schema, "khash", 1) ;
  if (!of) die ("failed to open syncmer file %s to read", syncSeqHashFileName) ;
  
  // read the syncmer hash parameters
  while (oneReadLine (of) && of->lineType != 'h') ;
  if (of->lineType == 'h') checkParams (of) ;
  else die ("sync file %s has no 'h' parameters record", syncSeqHashFileName) ;

  // read in the kmerhash
  khseq = kmerHashReadOneFile (of) ;
  if (khseq->len != params.w + params.k)
	  die ("syncmer len mismatch %d != %d + %d", khseq->len, params.w, params.k) ;

  // read in the node counts 
  if (of->lineType != 'C' || oneLen(of) != kmerHashMax(khseq))
    die ("failed to find expected C line in khash file") ;
  cseq = new (khseq->max, I64) ;
  memcpy (cseq, oneIntList(of), sizeof(I64) * khseq->max) ;

  fprintf (stderr, "read %llu nodes from %s\n", kmerHashMax(khseq), syncSeqHashFileName) ;
  oneFileClose (of) ; // close the khash file

  of = oneFileOpenRead (syncSeqFileName, schema, "syncseq", 1) ;
  if (!of) die ("failed to open syncmer file %s to read", syncSeqFileName) ;

  // read the syncmer hash parameters
  while (oneReadLine (of) && of->lineType != 'h') ;
  if (of->lineType == 'h') checkParams (of) ;
  else die ("sync file %s has no 'h' parameters record", syncSeqFileName) ;

#ifdef PRINT_SYNC_PROFILE
  I64 *rcnts = new (of->info['S']->given.max, I64);
  I64 *scnts = new (of->info['S']->given.max, I64);
  I64 *sposi = new (of->info['S']->given.max, I64);
#endif

  // go to the first sequence
  while (oneReadLine (of) && of->lineType != 'S') ;
  while (of->lineType == 'S')
    { I64 i, n, nr, ns, index, low = 0, high = 0, *syncs = oneIntList(of) ;
      bool exist;
      double c, w = .0 ;
      for (i = 0; i < oneLen(of); i++)
        { ns = cseq[syncs[i]-1] ;
          exist = kmerHashFindPacked (khref, packseq(khseq, syncs[i]), &index) ;
          nr = exist? cref[index-1] : 0 ;
          c = (double) nr / ns ;
          if (nr < minC)
            low += 1 ;
          else if (nr <= maxC)
            high += 1 ;
          w += c ;
#ifdef PRINT_SYNC_PROFILE
          scnts[i] = ns;
          rcnts[i] = nr;
#endif
        }
      n = oneLen(of) ;
#ifdef PRINT_SYNC_PROFILE
      while (oneReadLine (of) && of->lineType != 'P') ;
      if (of->lineType != 'P') die ("sync file %s has no expected 'P' line", syncSeqFileName) ;
      memcpy (sposi, oneIntList(of), sizeof(I64) * n); 
#endif
      while (oneReadLine (of) && of->lineType != 'R') ;
      if (of->lineType != 'R') die ("sync file %s has no expected 'R' line", syncSeqFileName) ;
#ifdef PRINT_SYNC_PROFILE
      for (i = 0; i < n; i++)
        fprintf(stderr,"PROF %9lld %9lld %12lld %9lld %9lld %9.1f\n",oneInt(of,1),i,sposi[i],scnts[i],rcnts[i],(double)rcnts[i]/scnts[i]);
#endif
      exist = (high >= minN && (double) low / n <= maxF && w / n >= minA) ;
      fprintf(out, "%-9lld %-9lld %-9lld %-9lld %9.3f %s\n", oneInt(of, 1), n, low, high, w / n, exist? "KEEP" : "DROP") ;
      if (exist) nSeq += 1 ;
      if (!oneReadLine (of)) break ;
    }
  oneFileClose (of) ; // close the khash file

  fprintf(stderr, "extracted %lld reads\n", nSeq) ;

#ifdef PRINT_SYNC_PROFILE
  free (rcnts);
  free (scnts);
  free (sposi);
#endif

  if (out != stdout)
    fclose(out) ;
  
  free (cref) ;
  free (cseq) ;

  kmerHashDestroy (khref) ;
  kmerHashDestroy (khseq) ;

  oneSchemaDestroy (schema);
	
  free(getCommandLine()) ;

  timeUpdate (stderr) ;
  fprintf (stderr, "total: ") ; timeTotal (stderr) ;
  return 0 ;
}

/*********************** end of file **********************/
