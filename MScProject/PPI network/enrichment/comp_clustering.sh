#!/bin/sh

echo "Running cluster enrichment comparsisons"

echo "\n--------\n"
echo "running BP vs Disease"
./run -opt 2 -Comfile ../clustering/Clustering/PPI_Network/sgG5_communities.csv -Annofile Annotations/flatfile.go.BP.csv -Annofile Annotations/flatfile_human_gene2HDO.csv -ext BPvsDisb -o OUT/Comparisions

echo "\n--------\n"
echo "running BP vs DE"
./run -opt 2 -Comfile ../clustering/Clustering/PPI_Network/sgG5_communities.csv -Annofile Annotations/flatfile.go.BP.csv -Annofile Annotations/flatfile_DE.csv -ext BPvsDE -o OUT/Comparisions

echo "\n--------\n"
echo "running Dis vs DE"
./run -opt 2 -Comfile ../clustering/Clustering/PPI_Network/sgG5_communities.csv -Annofile Annotations/flatfile_human_gene2HDO.csv -Annofile Annotations/flatfile_DE.csv -ext DEvsDis -o OUT/Comparisions

echo "\n--------\n"
echo "running BP vs Dis vs DE"
./run -opt 4 -Comfile ../clustering/Clustering/PPI_Network/sgG5_communities.csv -Annofile Annotations/flatfile.go.BP.csv -Annofile Annotations/flatfile_human_gene2HDO.csv -Annofile Annotations/flatfile_DE.csv -o OUT/Comparisions



