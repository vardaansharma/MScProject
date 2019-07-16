#!/bin/bash

#run enrichment for spectral clustering
echo running emrichment of spectral
echo human_gene2HDO
./run -Comfile ../clustering/Clustering/PPI_Network/spectral_communities.csv -Annofile Annotations/flatfile_human_gene2HDO.csv -opt 1 -setFDR BY -onesided -noPerm -o OUT/Spectral -ext human_gene2HDO

echo interpro family
./run -Comfile ../clustering/Clustering/PPI_Network/spectral_communities.csv -Annofile Annotations/flatfile.interpro.Family.csv -opt 1 -setFDR BY -onesided -noPerm -o OUT/Spectral -ext interpro_family

echo interpro domain
./run -Comfile ../clustering/Clustering/PPI_Network/spectral_communities.csv -Annofile Annotations/flatfile.interpro.Domain.csv -opt 1 -setFDR BY -onesided -noPerm -o OUT/Spectral -ext interpro_domain


#run enrichment for geo clustering
echo running enrichment of geo
echo human_gene2HDO
./run -Comfile ../clustering/CDMSuite_cpp/Geo/communityout.txt -Annofile Annotations/flatfile_human_gene2HDO.csv -opt 1 -setFDR BY -onesided -noPerm -o OUT/geo -ext human_gene2HDO

echo interpro_family
./run -Comfile ../clustering/CDMSuite_cpp/Geo/communityout.txt -Annofile Annotations/flatfile.interpro.Family.csv -opt 1 -setFDR BY -onesided -noPerm -o OUT/geo -ext interpro_family

echo interpro domain
./run -Comfile ../clustering/CDMSuite_cpp/Geo/communityout.txt -Annofile Annotations/flatfile.interpro.Domain.csv -opt 1 -setFDR BY -onesided -noPerm -o OUT/geo -ext interpro_domain


#run enrichment for sgG1 clustering
echo runnning enrichment of sgG1
echo human_gene2HDO
./run -Comfile ../clustering/Clustering/PPI_Network/sgG1_communities.csv -Annofile Annotations/flatfile_human_gene2HDO.csv -opt 1 -setFDR BY -onesided -noPerm -o OUT/sgG1 -ext human_gene2HDO

echo interpro_family
./run -Comfile ../clustering/Clustering/PPI_Network/sgG1_communities.csv -Annofile Annotations/flatfile.interpro.Family.csv -opt 1 -setFDR BY -onesided -noPerm -o OUT/sgG1 -ext interpro_family

echo interpro domain
./run -Comfile ../clustering/Clustering/PPI_Network/sgG1_communities.csv -Annofile Annotations/flatfile.interpro.Domain.csv -opt 1 -setFDR BY -onesided -noPerm -o OUT/sgG1 -ext interpro_domain

#run enrichment for lourvain clustering
echo running enrichment of lourvain
echo human_gene2HDO
./run -Comfile ../clustering/Clustering/PPI_Network/louvain_communities.csv -Annofile Annotations/flatfile_human_gene2HDO.csv -opt 1 -setFDR BY -onesided -noPerm -o OUT/louvain -ext human_gene2HDO

echo interpro_family
./run -Comfile ../clustering/Clustering/PPI_Network/louvain_communities.csv -Annofile Annotations/flatfile.interpro.Family.csv -opt 1 -setFDR BY -onesided -noPerm -o OUT/louvain -ext interpro_family

echo interpro domain
./run -Comfile ../clustering/Clustering/PPI_Network/louvain_communities.csv -Annofile Annotations/flatfile.interpro.Domain.csv -opt 1 -setFDR BY -onesided -noPerm -o OUT/louvain -ext interpro_domain

#run enrichment for fc clustering
echo running enrichment of fcc
echo human_gene2HDO
./run -Comfile ../clustering/Clustering/PPI_Network/fc_communities.csv -Annofile Annotations/flatfile_human_gene2HDO.csv -opt 1 -setFDR BY -onesided -noPerm -o OUT/fc -ext human_gene2HDO

echo interpro_family
./run -Comfile ../clustering/Clustering/PPI_Network/fc_communities.csv -Annofile Annotations/flatfile.interpro.Family.csv -opt 1 -setFDR BY -onesided -noPerm -o OUT/fc -ext interpro_family

echo interpro domain
./run -Comfile ../clustering/Clustering/PPI_Network/fc_communities.csv -Annofile Annotations/flatfile.interpro.Domain.csv -opt 1 -setFDR BY -onesided -noPerm -o OUT/fc -ext interpro_domain
