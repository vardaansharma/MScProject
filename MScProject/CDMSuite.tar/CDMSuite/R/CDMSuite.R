
CDMSuite <- function(
                     seed=1,
                     algorithm=1,
                     weighted="nw",
                     network){
  .Call( "CDMSuite",
        as.integer(seed),
        as.integer(algorithm),
        weighted,
        network,
        PACKAGE = "CDMSuite" )
}


help <- function(){

   print("CDMSuite requires 4 arguments: ");
   print("argument 1: random number seed  ");
   print("argument 2: the type of algorithm to run ");
   print("          : 1 = Geodesic edge Betweenness");
   print("          : 2 = Random edge Betweenness");
   print("          : 3 = Spectral Betweenness");
   print("argument 3: specify if network file is *weighted or not: ");
   print("          : w  = Using a weighted network file ");
   print("          : nw = Using a non-weighted network file ");
   print("argument 4: the network file to run  ");
   print("Example   : ./run 1 1 testData/Kirate.wpairs ");
   print("*         : The structure of the network file is:");
   print("w         : A (interactor) \t B (interactor) \t W (weight) ");
   print("nw        : A (interactor) \t B (interactor) ");
   print("          : Where A and B are integers and W is a double.");
  
}
