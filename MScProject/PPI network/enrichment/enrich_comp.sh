#!/bin/sh

echo "runnning the BP vs Dis min 5"
./run -opt 2 -Comfile ../clustering/Clustering/PPI_Network/sgG5_communities.csv -Annofile Annotations/flatfile.go.BP.csv -Annofile Annotations/flatfile_human_gene2HDO.csv -ext BPvsDisb -o OUT/Comparisions/min5 -setFDR BY -minOV 5 -minOV 5 -minOV 5 -maxSS -Chi2

./run -opt 2 -Comfile ../clustering/Clustering/PPI_Network/sgG5_communities.csv -Annofile Annotations/flatfile.go.BP.csv -Annofile Annotations/flatfile_human_gene2HDO.csv -ext BPvsDisb -o OUT/Comparisions/min3 -setFDR BY -minOV 3 -minOV 3 -minOV 3 -maxSS -Chi2


echo "running BP vs DE"
./run -opt 2 -Comfile ../clustering/Clustering/PPI_Network/sgG5_communities.csv -Annofile Annotations/flatfile.go.BP.csv -Annofile Annotations/flatfile_DE.csv -ext BPvsDE -o OUT/Comparisions/min3 -setFDR BY -minOV 3 -minOV 3 -minOV 3 -maxSS -Chi2

./run -opt 2 -Comfile ../clustering/Clustering/PPI_Network/sgG5_communities.csv -Annofile Annotations/flatfile.go.BP.csv -Annofile Annotations/flatfile_DE.csv -ext BPvsDE -o OUT/Comparisions/min5 -setFDR BY -minOV 5 -minOV 5 -minOV 5 -maxSS -Chi2

echo "running Dis vs DE"
./run -opt 2 -Comfile ../clustering/Clustering/PPI_Network/sgG5_communities.csv -Annofile Annotations/flatfile_human_gene2HDO.csv -Annofile Annotations/flatfile_DE.csv -ext DEvsDis -o OUT/Comparisions/min5 -setFDR BY -minOV 5 -minOV 5 -minOV 5 -maxSS -Chi2

./run -opt 2 -Comfile ../clustering/Clustering/PPI_Network/sgG5_communities.csv -Annofile Annotations/flatfile_human_gene2HDO.csv -Annofile Annotations/flatfile_DE.csv -ext DEvsDis -o OUT/Comparisions/min3 -setFDR BY -minOV 3 -minOV 3 -minOV 3 -maxSS -Chi2
