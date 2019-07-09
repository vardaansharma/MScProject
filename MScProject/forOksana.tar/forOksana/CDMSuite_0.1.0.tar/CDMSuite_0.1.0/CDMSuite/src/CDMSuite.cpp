#include <Rcpp.h>
#include "Headers.h"
#include "network.h"
#include "readfile.h"
#include "SpectralModularity.h"

// [[Rcpp::export]]
Rcpp::List spectral( Rcpp::DataFrame     DF,
		     Rcpp::IntegerVector CnMIN=1,
		     Rcpp::NumericVector TOL=0.00001 ){
  //Rcpp::List spectral( Rcpp::List params ){

  //For more information wrapping and packaging C/C++ in R see:
  //[1] https://www.gormanalysis.com/blog/exposing-a-cpp-student-class-with-rcpp/ (building R package and adding your C/C++ code)
  //[2] https://www.youtube.com/watch?v=DWkIbk_HE9o (setting up roxygen in R package)
  //[3] http://web.mit.edu/insong/www/pdf/rpackage_instructions.pdf
  //[4] http://r-pkgs.had.co.nz/src.html
  
  //Development steps
  //1) Open RStudio
  //2) edit/change code
  //3) Run pkgbuild::compile_dll() to compile your C++ code
  //4) Run devtools::document() to automatically build the NAMESPACE file package documentation.
  //5) Run devtools::load_all() to compile the C++ code and load our package

  //To build and install package into R
  //6) cd /afs/inf.ed.ac.uk/user/c/cmclean5/WORK/STUDIES
  //7) Run R CMD build mynewpackage
  //8) Run R CMD INSTALL mynewpackage_0.1.tar.gz
  //9) Start R
  //10) library(CDMSuite)
  //11) CDMSuite::spectral(...)
  
  int i,j,KK;

  int Cn_min       = 1;
  double tol       = 0.00001;
  int N            = 0;
  int M            = 0;
  double *A        = 0;
  edgelist *el     = 0;
  bool useLoops    = false;
  bool checkM      = true;
  
  //initialise network
  readfile *reader          = 0;
  network *gg               = new network();
  SpectralModularity *model = 0;
  
  
  int ncols = DF.length();
  int nrows = DF.nrows();

  if( (ncols > 0) && (nrows > 0) ){

    if( CnMIN.length() == 1 ){
      if( (CnMIN[0] > 0) ){
	Cn_min = CnMIN[0];
      }
    }
      
    if( (TOL.length() == 1) ){
      if( (TOL[0] > 0) ){
	  tol = TOL[0];
	}
    }  
    
    KK              = nrows*ncols;
    string *DATASET = new string[KK];
      
    for( i=0; i<ncols; i++ ){
      Rcpp::StringVector col = DF[i];
      for( j=0; j<nrows; j++ ){
	DATASET[(j*ncols)+i] = col[j];
      }
    }
  
    //load edgelist into network    
    reader = new readfile( gg, DATASET, ncols, nrows );

    //build Adjaceny Matrix
    gg->buildNetworkReps( useLoops, checkM );
    //---
      
    N = gg->getN();
    M = gg->getM2();
    A = gg->getA();

    if( N != 0 && M != 0 ){

      //set-up clustering alg.
      model = new SpectralModularity(gg,el,A,N,M);	
      model->settol( tol );
      model->setMinCn( Cn_min );
	
      //--- run spectral clustering
      int cal = model->calculateSpectralModularity();
	
      //--- reorder community numbers in network
      gg->reorderK();        

    }  

  }  
  
  
  if( gg->getN() > 0 ){

    N = gg->getN();
    
    //--- output the node label and its cluster
    Rcpp::StringVector  ID   (N);
    Rcpp::NumericVector Coms (N);
    
    for( i=0; i<N; i++ ){
      ID[i]   = gg->V[i].label;
      Coms[i] = gg->V[i].K;
    }

    // Create a named list with the above quantities
    return Rcpp::List::create(Rcpp::Named("ID") = ID,
			      Rcpp::Named("K")  = Coms);
      
  } else {

    //--- output the node label and its cluster
    Rcpp::StringVector  ID   (0);
    Rcpp::NumericVector Coms (0);
    
    // Create a named list with the above quantities
    return Rcpp::List::create(Rcpp::Named("ID") = ID,
			      Rcpp::Named("K")  = Coms);
    
  } 

  
}//spectral
  
