//#if !defined(NETWORKENRICHMENT_INCLUDED)
//#define NETWORKENRICHMENT_INCLUDED
#include "Headers.h"
#include "buildSets.h"

//define guards, so headers are declare only once.
#ifndef NETWORKENRICHMENT_H
#define NETWORKENRICHMENT_H

class NetworkEnrichment : buildSets {


 public:
   NetworkEnrichment();
   NetworkEnrichment( vector<string> ); 
  ~NetworkEnrichment();

  int   getKOffset();
  void  setKOffset( int = 1 );
  
  double getMINOVERLAP ( int );
  void   setMINOVERLAP ( int , int  );
  void   setANNOindex( int );
  void   setFDRmethod ( const char* );
  int    calculateOverlapinCommunities(bool,          const char*, const char*, bool = false, bool = false );
  int    calculateOverlapinCommunities(int, int,      const char*, const char*, bool = true );
  int    calculateOverlapinNetwork    (int, int,      const char*, const char*, bool = true );
  int    calculateOverlapinNetwork    (int, int, int, const char*, const char*, bool = true );
  int    calculateOverlapinCommunities(int,           const char*, const char*, bool = true);
  
 private:
  void               freeMemory();
  void               setSeed(bool=false, int=0);
  unsigned long int  getSeed();

  void geneAssociations( int, int, int, int &, int & );
  void geneAssociations( int, int, int[] );
  void overlapinNetwork();
  void overlapinNetork( int, int );
  void overlapinNetork( int, int, int );
  void permutation( double );
  double prob_overlap( int, int, int, int );//Hypergeometric distribution
  double prob_overlap( int, int, int, int, int );//intersection between three sets

  void overlapinComsHypergeometricTestRnd();
  void overlapinComsHypergeometricTest   ();
  void overlapinComsHypergeometricTest   (int, int);
  void overlapinNetHypergeometricTest    (int, int);
  void overlapinNetHypergeometricTest    (int, int, int);
  bool _min( double, double );
  bool _max( double, double );
  void CalculateFDR_BY( vector<pairDoubInt>, double, double &, double &, int &, vector<pairDoubInt> &);
  void CalculateFDR_BH( vector<pairDoubInt>, double, double &, double &, int &, vector<pairDoubInt> &);
  void CalculateFDR_BL( vector<pairDoubInt>, double, double &, double &, int &, vector<pairDoubInt> &);
  void calculateFDR( int=1, int=-1, int=-1, int=-1 );
  
  void printOverlapinCommunities   (const char *, const char *, bool);
  void printOverlapinCommunitiesAlt(const char *, const char *, bool);
  void printOverlapinCommunities   (int, int, const char *, const char *, bool);
  void printOverlapinNetwork       (int, int, const char *, const char *, bool);
  void printOverlapinNetwork       (int, int, int, const char *, const char *, bool);
  void printFDR                    ();
  
  //GSL random number and seed
  unsigned long int seed;
  gsl_rng *g;
  
  //set output file delineator(s), in print functions
  static const int DELSIZE = 1;
  char dels[DELSIZE]; //use for output file delineator(s
  
  vector<string> INfiles;  //store input file paths.
  vector<string> baseNAME; //base names for print function, can change.
  //vector<string> setNames;
  vector<string> FDRmethods;
  
  vector<bool> freedMemory;//flags to indicate if memory needs freeing 
  int ANNOindex;           //Selected annotation set index
  int NoP;                 //No: of permutation studies to run  
  int FDRtest;             //Select which FDR test to run
  int KOFFSET;             //Value added to community number  
  bool isOFFSET;           //Flag if we have offset community numbers

  double FDR;
  double PV;
  double SIGMA;
  int    LEVEL;
  int    TESTS;
  
  
  //global variable values
  static const int MAXEXPO     = 700; //see prob_overlap functions
  static const int SIGMASIZE   = 3;
  static const int OVERLAPSIZE = 4; 
  static const int BUFFERSIZE  = 250;
  
  double sigma[SIGMASIZE];
  double bonferroni[SIGMASIZE];
  double MINOVERLAP[OVERLAPSIZE];

  //arrays for Hypergeometric tests  
  int* comSIZE;
  int* geneCOM;
  int* annoSIZE;
  
  double* studies;
  double* p_values;
  double* permute;
  double* overlap;
  double* muab;

  double* padjusted;

};

#endif

