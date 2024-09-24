1. add checks if requested hmm exists - done
2. add the possibility to run on a provided fasta file - e.g. from broccoli - done
3. phylogeny module - so far local - done
4. make the search outputs obvious from the command line name

python main.py search -f data/sample.fasta -g data/genefam.tsv Myosin # creates results_annotation/searches/myo.Myosin.domains.fasta

mkdir -p results_annotation/alignments
python main.py align -f results_annotation/searches/myo.Myosin.domains.fasta -o results_annotation/alignments/test.aln -c 10
mkdir -p results_annotation/gene_trees
python main.py phylogeny -f results_annotation/alignments/test.aln -o results_annotation/gene_trees/test -c 15

