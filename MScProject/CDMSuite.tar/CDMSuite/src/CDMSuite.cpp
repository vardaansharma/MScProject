
#include <cmath>
#include <Rcpp.h>
#include "files/Headers.h"
#include "files/communityDetection.cpp"
#include <Rdefines.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

using namespace Rcpp ;

extern "C"{
  SEXP CDMSuite( SEXP _seed, 
		 SEXP _a_type,
		 SEXP _is_weighted,
		 SEXP _networkfile
		 ){

    //--- The R Rcpp object to return 
    List ret; 	

    //--- Convert input variables from R to C
    CharacterVector _w_(_is_weighted);
    CharacterVector _f_(_networkfile);

    int  _s        = INTEGER(_seed)[0];
    int  _a        = INTEGER(_a_type)[0];    
    const char *_w = _w_[0];
    const char *_f = _f_[0];
    string __w     = _w;

    //--- Check input variables

    //--- Check the algorithm to use
    if( _a < 1 || _a > 3 ){
      cout << "argument 2: the type of algorithm to run needs to be either (1,2,3): " << endl;
      cout << "          : 1 = Geodesic edge Betweenness" << endl;
      cout << "          : 2 = Random edge Betweenness"   << endl;
      cout << "          : 3 = Spectral Betweenness"      << endl;
      return(ret);
  }

    //--- Check if using a weighted or unweighted network.
    if( __w.compare("w") == 0 ){
      ;
    } else {
      if( __w.compare("nw") == 0 ){
	;
      } else {
	cout << "argument 3: specify if network file is weighted or not: " << endl;
      cout << "          : w  = Using a weighted network file " << endl;
      cout << "          : nw = Using a non-weighted network file " << endl;    
      return(ret);
      }
    }

    //--- Check input file is valid
    ifstream fscan(_f,ios::in);
    if( fscan.is_open() != true ){
      cout << "> Error opening " << _f << endl;
      cout << "> Use files absolute path name " << endl;
      return(ret);
    }

    //--- Execute the community detection algorithm
    communityDetection *cdm = new communityDetection( _s, _a, _w, _f );


    //--- Return C objects for consensus matrix
    vector<int> keyi = cdm->getKey_listi();
    vector<int> keyj = cdm->getKey_listj();
    vector<int> keyk = cdm->getKey_listk();


    //--- Convert C consensus object into R
    List consensus = List::create(Named("key_listi") = keyi,
				  Named("key_listj") = keyj,
				  Named("key_listk") = keyk);


    //--- Return C objects for each community and it's size.
    vector<int> com_list = cdm->getComList();
    vector<int> com_size = cdm->getComsize();

    
    //--- Convert C community objects into R
    List comout    = List::create(Named("No_communities") = com_list,
				  Named("Community_size")  = com_size); 


    //--- Return C objects for Network info
    int            node_size = cdm->getNodeSize();
    int            edge_size = cdm->getEdgeSize();    
    double         cpu_time  = cdm->getCPUtime();    


    //--- Convert C Network info objects into R
    List NtwkInfo; 
    NtwkInfo = List::create(Named("No_Nodes")    = node_size,
			    Named("No_Edges")    = edge_size,
			    Named("CPU_time")    = cpu_time,
			    Named("Communities") = comout);


    //--- C objects for the Modularity
    vector<double> mod;
    vector<double> mod_err;    
    
    double         mod_max;
    double         mod_err_max;
    double         mod_spec;


    List EdgeBet;
    List Spec;
    if( _a != 3 ){

      //--- Return C Modularity objects for Geodesic or RandomWalk 
      //--- Edge Betweenness.
      mod         = cdm->getMod();
      mod_err     = cdm->getModErr();

      mod_max     = cdm->getModMax();
      mod_err_max = cdm->getModErrMax();

      //--- Convert C Modularity objects for Geodesic or RandomWalk 
      //--- into R
      EdgeBet = List::create(Named("Mod_max")          = mod_max,
			     Named("Mod_err_max")      = mod_err_max,
			     Named("Modularity")       = mod,
			     Named("Modularity_error") = mod_err);

      //--- Store in R results object
      ret = List::create(Named("NetworkInfo") = NtwkInfo,
			 Named("Modularity")  = EdgeBet,
			 Named("consensus")   = consensus );    
     
    } else {

      //--- Return C Modularity objects for Spectral algorithm
      mod_spec    = cdm->getModSpectal();

      //--- Convert C Modularity objects for Spectral algorithm into R
      Spec = List::create(Named("Spectral_Modularity") = mod_spec);

      //--- Store in R results object
      ret = List::create(Named("NetworkInfo") = NtwkInfo,
			 Named("Modularity")  = Spec,
			 Named("consensus")   = consensus );    

    }

    //--- Remove C community detection object 
    cdm->~communityDetection();


    //--- Return R results object
    return(ret);

  }

}//extern C
