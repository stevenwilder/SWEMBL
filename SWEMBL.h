#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <ctype.h>
#include <zlib.h>
#include <time.h>

struct region
{ 
  char chr[50];
  long int start, end, seen_start, seen_end;
  unsigned int pos;
  double count, refcount, total, print_total, summit, peakheight, peakheightmax;
};

struct startendcount
{
  long int startend[2];
  double count;
};

struct countend
{
  double count;
  long int end;
};

struct readinfo
{
  char chr[50], prevchr[50];
  long int pos, prevpos, negpos;
  char strand[2];
  int qual, pairlength, pairflag, changechr;
  double count, pvecount, negcount, refcount;
};

typedef enum {false, true}  boolean;

struct elem {
  long int d[2];
  double s;
  struct elem *next;
};

typedef struct elem         elem;

struct stack {

  int         cnt;
  elem        *head;
  elem        *tail;
};

typedef struct stack        stack;

struct regionelem {
  struct region     r;
  struct regionelem *next;
};

typedef struct regionelem         regionelem;

struct regionstack {
  int         cnt;
  regionelem  *head;
  regionelem  *tail;
};

struct chrelem{
  char     chr[50];
  struct   stack *stk;
  struct chrelem *next;
};

typedef struct chrelem   chrelem;

struct chrstack {
  int      cnt;
  chrelem  *head;
  chrelem  *tail;
};

typedef struct chrstack  chrstack;

typedef struct regionstack    regionstack;

struct param {
  double bg ,posbg, longbg, qualcutoff, resultcutoff, refqualcutoff, refpen, maxtotal, proportion, threshold, refthreshold, binom_p;
  unsigned int min_above_bg, seqlength, pen_inc, usequal, refseqlength, killifnoref, totalreads, totalrefreads;
  int fraglength, reffraglength, paired, ref, wig, contigrows, comp, quietmode, overlapmode, bootstrap, newpeakversion;
  long int nlines;
};

typedef struct param   param;


/*struct savedreflengths {

   char                 chr[50];
   unsigned long int    pvecount;
   unsigned long int    negcount;

   } ; */

// typedef struct savedreflengths savedreflengths;

struct contigconv {

  char contig[50];
  char chr_random[50];
  long int start,end;
};

typedef struct contigconv   contigconv;

char filetype[2];
char filemode[2];
char pvestring[2];
char negstring[2];
int endcol;
int zip;

struct comparisons {
  int print, mode, nseqregions, nseqregionsoverlap, ncompregions, ncompregionsoverlap;
  long int lenseqregions, lenoverlap, lencompregions, compprioroverlap;
};


typedef struct comparisons comparisons;

struct prevregions {
  struct region last_but_1_region;
  char compchr[50];
  long int compstart, compend;
};

typedef struct prevregions prevregions;

FILE *in_fp;
FILE *out_fp;
FILE *ref_fp;
FILE *wig_fp;
FILE *NT_fp;
FILE *comp_fp;
FILE *overcomp_fp;
FILE *overreg_fp;

struct comparisons overlap;
struct prevregions prev;

unsigned int rseed;
double ppois[10];
double refppois[10];

//int main(int argc, char **argv);
void store_results(struct readinfo read, struct region possend, struct param par, struct region *lastregion);
struct param get_parameters(char infile[1000], char outfile[1000], char reffile[1000], char wigfile[1000], char contigfile[1000], char compfile[1000], char compoverlapfile[1000], char regoverlapfile[1000], int argc, char **argv);
void print_help();
void open_files(int ref, int wig, int comp, char *infile, char *outfile, char *reffile, char *wigfile, char *contigfile, char *compfile, char *compoverlapfile, char *regoverlapfile);
int platform_strand_endline();
struct param get_firstline(int *endfile, char *nextline, struct param par);
char *egets(char *line, int n, FILE *fp);
void print_output_header(struct param par, char *infile, char *reffile);
void readrefline(struct readinfo *ref, int *refendfile, int reffraglength, int paired, int bootstrap);
void readsampleline(struct readinfo *read, char *nextline, int *endfile, int fraglength, int seqlength, int qualcutoff, int paired, int overlapmode, int bootstrap) ;
long int set_negpos(struct readinfo *read, int *endfile, int seqlength, int fraglength);
void set_fragpos(struct readinfo read, long int *fragpos, long int *fragend, long int negpos, long int fraglength);
void reference_reads(struct readinfo read, struct readinfo *refread, int *refendfile, struct param par, int *nsavedrefchrs, long int *refnegpos, long int negpos, long int *fragpos);
void end_chr(struct readinfo read, struct region *curr, struct param par, int fragpos, int fragend, struct region *lastregion);
void ref_endchr(struct readinfo read, struct region *curr, struct param par, int fragpos, int fragend, struct region *lastregion);
void update_total(int dir, struct region *curr, struct region *possend, long int negpos, long int fragpos, long int fragend, long int prevfragpos, struct param par, struct readinfo read, struct region *lastregion);
void update_fragpos(long int *fragpos, long int seqlength, long int fraglength, struct readinfo read);
void update_fragpos_ref(long int *fragpos, char readchr[], char refreadchr[]);
int separate(char strnd[2], int pairlen, long int pos, int bootstrap);
struct region new_region(struct readinfo read, struct param par, long int fragpos, long int fragend, int dir);
int set_back_fragpos(long int *fragpos, long int *fragstart);
double sample_back_count(long int fragpos, long int *fragstart, double threshold, stack *stk);
double ref_back_count(long int fragpos, double refthreshold, stack *stk);
void output_regions(struct region possend, struct region *olr, struct region pvestrandregion, struct param par);
struct region check_and_print_regions(struct region olr, struct region *lastregion, int comp, double binom_p, int peakversion);
void print_out_lastregion(struct region *lastregion, double binom_p, int peakversion);
void establish_wiggle(char infile[1000], char contigfile[1000], struct contigconv conversion[], int *wigonevalue);
void push_wigcounts(struct readinfo read, long int fragpos, long int prevfragpos,long int fragend,int fraglength, struct contigconv conversion[], int contigrows, char wigchr[], long int *wigstart, int paired, int *wigonevalue);
void export_wiggle(contigconv conv[], int n, char wigchr[], long int *wigstart, int fraglength, int paired, int wigonevalue);
//int wiggle_paired_one_value(stack *stk);
//void store_comp_locations(struct stack *compresults[], char *compchrs[], int *compfound[]);
void compare_lastregion(struct region *lastregion);
int store_comp_location();
int seenchr(char chr[], struct chrstack *cstk);
void print_out_overlap();
void nearestregion(struct chrstack *cstk, struct region *lastregion);

void        add_to_position(double n, int posn, stack *stk);
void add_to_or_insert_start_position(long int start, long int end, double n, stack *stk);
long int return_position(int posn, stack *stk);
void        initialize(stack *stk);
void        push(long int start, long int end, double count, stack *stk);
void        unshift(long int start, long int end, double count, stack *stk);    //push in book
struct countend cntequalhead(long int pos, stack *stk);
void sortunshift(long int start, long int end, double count, stack *stk);
void assign(stack *stk, stack *istk);
double median(stack *stk);
int cnt(stack *stk);
struct startendcount        shift(stack *stk);        //pop in book
struct startendcount        head(stack *stk);        //top in book
struct startendcount        tail(stack *stk);       //not in book
boolean     empty(const stack *stk);
void        initializeregion(regionstack *stk);
int regioncnt(regionstack *stk);
void        unshiftregion(struct region r, regionstack *stk);    //push in book
struct region   shiftregion(regionstack *stk);        //pop in book
struct region   headregion(regionstack *stk);        //top in book
boolean     emptyregion(const regionstack *stk);

void initializechr(chrstack *stk);
int chrcnt(chrstack *stk);
void pushchr(long int start, long int end, double count, chrstack *cstk);
void addchr(char chr[50], chrstack *cstk);
struct startendcount shiftchr(chrstack *cstk);
struct stack *removechr(chrstack *chrstk);
struct startendcount headchr(chrstack *cstk);
struct startendcount nextheadchr(chrstack *cstk);
struct stack *firstchr(chrstack *stk);
struct chrelem *firstchrelem(chrstack *cstk);
struct chrelem *secondchrelem(chrstack *cstk);
struct stack *lastchr(chrstack *cstk);
struct chrelem *lastchrelem(chrstack *cstk);
boolean emptychr(const chrstack *stk);
long int count_sample_lines(struct param *par, char infile[1000], int samplecount);
void push_peakheight(int dir, long int start, long int end, double score, struct region *curr, struct region *possend);
void calculate_peak_summit(struct region *curr, long int begin, long int finish, struct region *possend);

struct stack pvepos, backpvepos, backnvepos, refpvepos, refnvepos, savedpverefs, savednegrefs, savedrefpos, refbackpos, wigcounts, peakheights, peaksummits;
struct regionstack outregions;
struct chrstack refchrs_pve, refchrs_nve, compchrs;

void poiscalc(double ratio, double ppois[10], double refppois[10]);
int rpois(int bootstrap, int ref);
double pbinom(int q, int n, double p);
double pnorm(double x);
