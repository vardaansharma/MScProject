#ifndef COMMUNITYDETECTION_H
#define COMMUNITYDETECTION_H

#include "Headers.h"
#include "edge.h"
#include "node.h"
#include "Helper.h"
#include "readInputFile.h"
#include "nr3.h"
#include "ran.h"
#include "ludcmp.h"
#include "eigen_sym.h"

class communityDetection{

public:	

  communityDetection( int, int, const char*, const char*);
  ~communityDetection();

//--- METHOD DECLARATIONS
//-------------------------------------------------------------------------------------

//--- Methods for Geodesic edge Betweenness algorithm
//-------------------------------------------------------------------------------------
//--- Finding and evaluating community structure in networks, M. Newman and M. Girvan
//--- Physical Review E69, 026113 (2004). --- [1]
  void calculateEdgeBetweennessGeodesic();
  void assignNodeWeights( int );

//--- Methods for Random edge Betweenness algorithm
//-------------------------------------------------------------------------------------
//--- Finding and evaluating community structure in networks, M. Newman and M. Girvan
//--- Physical Review E69, 026113 (2004). --- [1]
 void calculateEdgeBetweennessRandom();
 void calculateRandomWalk( int, vector<node> );
 void setupMatrices  ();
 void upDateMatrices ();
 void getSubMatrix   ( int, vector<node> & );
 void removeMatrixRow( MatDoub Unr, MatDoub &outnr );
 void removeMatrixRow( MatDoub &out );
 void removeMatrixRow(int row, MatDoub &out);
 void addMatrixRows(MatDoub U, int c, MatDoub &Ui);
 void addMatrixRow(MatDoub U, int row, MatDoub &out);
 void visitedNodes( int );

//--- Methods for Spectral Modularity algorithm
//-------------------------------------------------------------------------------------
//--- Finding community structure in networks using eigenvectors of matrices
//--- M. Newman, arXix:physics/0605087v3 (2006). --- [2]
 void calculateSpectralModularity();
 void calculateB(MatDoub, MatDoub &);
 void calculateEigenVectors();
 void findLeadingEigenVectors(int & );
 void maximiseIndexVectors(int);
 void modifySplit(double, int );
 void maxModularity( double & );
 void deltaModularity( double &);
 void deltaModularityMax( int, double &);
 int  maxCommunity();
 
 void splitP( MatDoub, VecInt, int, double );
 void splitN( MatDoub, VecInt, int, double );

//--- Methods shared between algorithms
 void Modularity       ( double &, double & );
 int  findMaxModularity( vector<double> );
 int  deltaCommunity   ( int, int );
 int  delta            ( int, int );

 void getSubSample( vector<int> & );
 vector<node> storeNodes();
 int findMax( vector<double> );

 void printHelpMessage( const char *);
 void DrawProgressBar(int, double);

 //--- Methods to get output data
 vector<int>    getKey_listi();
 vector<int>    getKey_listj();
 vector<int>    getKey_listk();

 vector<double> getMod();
 vector<double> getModErr();
 
 double         getModMax();
 double         getModErrMax();
 double         getModSpectal();
 
 double         getCPUtime();
 int            getNodeSize();
 int            getEdgeSize();

 vector<int>    getComList();
 vector<int>    getComsize();

 void           printCommunityout();
 void           printConsensusout();

 private:

  //--- GLOBAL VARIABLES
  vector<node>  n;
  vector<node>  Gn;
  
  vector<edge>  elist;
  vector<edge>  Gelist;
  vector<edge>  totallist;
  
  vector<node>  emptyn;
  vector<edge>  emptye;

  vector<int> key_listi;
  vector<int> key_listj;
  vector<int> key_listk;

  int        *dataI;
  int        *datainI;
  double     *data;
  double     *datain;
  double     *datainO;
  double     *temp_score;
  stack<int> node_dist;

  readInputFile reader;
  Helper        ihelper;

  fstream* modularityscore;
  fstream* communityout;
  fstream* removededges;
  fstream* consensusout;
  
  vector<double>           vec_mod;
  vector<double>           vec_mod_err;
  vector<int>              vec_com_max;
  vector< vector<node> >   vec_nodes;

  vector<int> comlist;
  vector<int> comsize;

  MatDoub R;
  MatDoub Ri;
  MatDoub Rh;
  MatDoub A;
  MatDoub Ai;
  MatDoub Bi;

  VecInt C;
 
  MatDoub S;
  MatDoub V;
  MatDoub T;
  MatDoub Ti;
  MatDoub Rc;
  MatDoub Vi;

  MatDoub B;
  MatDoub Bm;

  MatDoub Bgi;
  
  VecInt keys_p;
  VecInt keys_n;

  MatDoub SI;
  VecDoub si;
  VecDoub visited;

  MatDoub u;
  VecDoub betai;

  double specQ; //Recorded Modularity value from the Spectral method
  double two_m; //2*m = Sum_i { K_i } = Sum_ij { A_ij }
  double _norm; //1.0/(4*m)
 
  int    com_max;

  double cpu_time_used;
  clock_t cstart; 
  clock_t cend;

  Ran _rand;
   
  int seed;
  int a_type;
  int w_type;
  string title;
  string if_weighted;
  string if_help;
  const char *file_network;
  const char *file_names;  


};

#endif
